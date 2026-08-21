# title: "ARID5B figure"
# author: "Sueanne Chear"


# siARID5B vs NTC analysis

library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(openxlsx)
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(grid)
library(pheatmap)
library(ggplotify)
library(clusterProfiler)
library(msigdbr)
library(ComplexUpset)
library(decoupleR)
library(OmnipathR)


#### Load counts + build sample_meta####

data     <- read.table("ARID5B_matrix.tsv", header = TRUE, sep = "\t")
metadata <- read.csv("ARID5B_metadata.csv")

count_cols <- grep("_count$", colnames(data), value = TRUE)
count_numbers <- sub("X5XV37Q_([0-9]+)_count", "\\1", count_cols)
count_numbers <- as.numeric(count_numbers)
new_names <- metadata$Sample_name[match(count_numbers, metadata$Number)]

counts <- data[, count_cols]
rownames(counts) <- data$gene_id
colnames(counts) <- new_names

# Map ENSEMBL -> SYMBOL and drop genes with no valid symbol
symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = rownames(counts),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

valid_symbol_genes <- names(symbol_map)[!is.na(symbol_map) & symbol_map != ""]
counts <- counts[valid_symbol_genes, ]
symbol_map <- symbol_map[valid_symbol_genes]

sample_meta <- data.frame(
  Sample = colnames(counts)
) %>%
  separate(
    Sample,
    into = c("Group", "Treatment", "Timepoint", "Rep", "Number"),
    sep = "_",
    remove = FALSE
  ) %>%
  unite(Replicate, Rep, Number, sep = "_") %>%
  mutate(
    Group     = factor(Group, levels = c("NTC", "siARID5B")),
    Treatment = factor(Treatment, levels = c("Basal", "E2", "MPA", "cAMP")),
    Timepoint = factor(Timepoint, levels = c("72h", "96h", "120h", "192h")),
    Condition = factor(
      paste(Group, Treatment, sep = "_"),
      levels = c(
        "NTC_Basal", "NTC_E2", "NTC_MPA", "NTC_cAMP",
        "siARID5B_Basal", "siARID5B_E2", "siARID5B_MPA", "siARID5B_cAMP"
      )
    ),
    FullCondition = factor(paste(Group, Treatment, Timepoint, sep = "_"))
  )

sample_meta <- sample_meta[match(colnames(counts), sample_meta$Sample), ]

treatments <- levels(sample_meta$Treatment)
timepoints <- levels(sample_meta$Timepoint)


#### DEG pipeline (filter -> normalize -> voom -> fit -> contrasts)####

dge <- DGEList(counts = counts)

design <- model.matrix(~0 + FullCondition, data = sample_meta)
colnames(design) <- levels(sample_meta$FullCondition)

keep <- filterByExpr(dge, design = design)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

v <- voom(dge, design, plot = FALSE)

####Save logCPM matrix + sample_meta####

fit <- lmFit(v, design)

contrast_df <- expand.grid(
  Treatment = treatments,
  Timepoint = timepoints,
  stringsAsFactors = FALSE
)

contrast_strings <- paste0(
  "siARID5B_", contrast_df$Treatment, "_", contrast_df$Timepoint,
  " - ",
  "NTC_", contrast_df$Treatment, "_", contrast_df$Timepoint
)

contrast_names <- paste0(
  "siARID5B_vs_NTC_", contrast_df$Treatment, "_", contrast_df$Timepoint
)

contr <- makeContrasts(contrasts = contrast_strings, levels = design)
colnames(contr) <- contrast_names

fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)


# Per-contrast OR + total-count-gate detectability filter

min_cpm          <- 1
min_total_count  <- 20

cpm_mat        <- cpm(dge)
raw_counts_mat <- dge$counts

all_ensembl <- rownames(fit2)
symbol_map_fit2 <- symbol_map[all_ensembl]

wb <- createWorkbook()

for (i in seq_along(colnames(contr))) {
  
  coef_name   <- colnames(contr)[i]
  treatment_i <- contrast_df$Treatment[i]
  timepoint_i <- contrast_df$Timepoint[i]
  
  ntc_samples <- sample_meta %>%
    filter(Treatment == treatment_i, Timepoint == timepoint_i, Group == "NTC") %>%
    pull(Sample)
  kd_samples <- sample_meta %>%
    filter(Treatment == treatment_i, Timepoint == timepoint_i, Group == "siARID5B") %>%
    pull(Sample)
  
  ntc_detected_both <- rowSums(cpm_mat[, ntc_samples, drop = FALSE] >= min_cpm) == length(ntc_samples)
  kd_detected_both  <- rowSums(cpm_mat[, kd_samples,  drop = FALSE] >= min_cpm) == length(kd_samples)
  
  ntc_passes_gate <- rowSums(raw_counts_mat[, ntc_samples, drop = FALSE]) >= min_total_count
  kd_passes_gate  <- rowSums(raw_counts_mat[, kd_samples,  drop = FALSE]) >= min_total_count
  
  keep_contrast <- (ntc_detected_both & ntc_passes_gate) | (kd_detected_both & kd_passes_gate)
  reliably_detected_ensembl <- rownames(cpm_mat)[keep_contrast]
  
  res <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "BH", sort.by = "P")
  res$ENSEMBL <- rownames(res)
  res <- res[res$ENSEMBL %in% reliably_detected_ensembl, ]
  
  res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")
  res$SYMBOL <- unname(symbol_map_fit2[res$ENSEMBL])
  res <- res[, c("SYMBOL", "ENSEMBL", setdiff(colnames(res), c("SYMBOL", "ENSEMBL")))]
  
  sheet_name <- substr(coef_name, 1, 31)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, res)
}

excel_file <- "siARID5B_vs_NTC_all_contrasts_limma_voom_filtered.xlsx"
saveWorkbook(wb, file = excel_file, overwrite = TRUE)


# Build all_filtered_res ONCE

contrast_sheets <- readxl::excel_sheets(excel_file)

parse_contrast <- function(sheet_name) {
  parts <- strsplit(sheet_name, "_")[[1]]
  n <- length(parts)
  list(Treatment = parts[n - 1], Timepoint = parts[n])
}

all_filtered_res <- map_dfr(contrast_sheets, function(sheet_name) {
  parsed <- parse_contrast(sheet_name)
  readxl::read_excel(excel_file, sheet = sheet_name) %>%
    mutate(
      Contrast  = sheet_name,
      Treatment = parsed$Treatment,
      Timepoint = parsed$Timepoint
    )
}) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatments),
    Timepoint = factor(Timepoint, levels = timepoints)
  )


####PCA + trajectory + sample correlation heatmap (uses v$E / sample_meta directly)####

logCPM <- v$E

pca <- prcomp(t(logCPM), scale. = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4],
  sample_meta
)

centroids <- pca_df %>%
  group_by(Treatment, Group, Timepoint) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop") %>%
  mutate(LineGroup = paste(Treatment, Group, sep = "_")) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("72h", "96h", "120h", "192h"), ordered = TRUE)) %>%
  arrange(LineGroup, Timepoint)

trajectory_segments <- centroids %>%
  group_by(LineGroup) %>%
  arrange(Timepoint, .by_group = TRUE) %>%
  mutate(PC1_next = lead(PC1), PC2_next = lead(PC2)) %>%
  filter(!is.na(PC1_next)) %>%
  ungroup()

arrow_segments <- trajectory_segments %>%
  mutate(
    x_start = PC1 + 0.38 * (PC1_next - PC1),
    y_start = PC2 + 0.38 * (PC2_next - PC2),
    x_end   = PC1 + 0.62 * (PC1_next - PC1),
    y_end   = PC2 + 0.62 * (PC2_next - PC2)
  )

treatment_colours <- c("Basal" = "#D55E00", "E2" = "#0072B2", "MPA" = "#009E73", "cAMP" = "#CC79A7")
timepoint_shapes  <- c("72h" = 16, "96h" = 17, "120h" = 15, "192h" = 18)
group_linetypes   <- c("NTC" = "solid", "siARID5B" = "dashed")

p_pca <- ggplot() +
  geom_point(data = pca_df, aes(x = PC1, y = PC2, colour = Treatment),
             size = 2, alpha = 0.45, show.legend = FALSE) +
  geom_path(data = centroids, aes(x = PC1, y = PC2, colour = Treatment, linetype = Group, group = LineGroup),
            linewidth = 0.9, alpha = 0.95, lineend = "round") +
  geom_segment(data = arrow_segments,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, colour = Treatment, linetype = Group),
               linewidth = 0.9, lineend = "round",
               arrow = arrow(length = unit(1.8, "mm"), type = "closed"), show.legend = FALSE) +
  geom_point(data = centroids, aes(x = PC1, y = PC2, colour = Treatment, shape = Timepoint),
             size = 4.2, stroke = 0.8) +
  scale_colour_manual(values = treatment_colours) +
  scale_shape_manual(values = timepoint_shapes, drop = FALSE) +
  scale_linetype_manual(values = group_linetypes) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    colour = "Treatment", shape = "Timepoint", linetype = "Group"
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 4, linewidth = 1)),
    shape  = guide_legend(order = 2, override.aes = list(colour = "black", size = 4)),
    linetype = guide_legend(order = 3, override.aes = list(colour = "black", linewidth = 1))
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.35),
    axis.title = element_text(size = 13, colour = "black"),
    axis.text  = element_text(size = 9, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    legend.key.width  = unit(1.1, "cm"),
    legend.key.height = unit(0.55, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

ggsave("PCA_centroid_trajectories.pdf", p_pca, width = 7.2, height = 4.8, units = "in", device = cairo_pdf)

group_colours <- c("NTC" = "#0072B2", "siARID5B" = "#D55E00")

p_pc34 <- ggplot(pca_df, aes(x = PC3, y = PC4, colour = Group)) +
  geom_point(size = 3.8, alpha = 0.80) +
  scale_colour_manual(values = group_colours) +
  labs(
    x = paste0("PC3: ", percentVar[3], "% variance"),
    y = paste0("PC4: ", percentVar[4], "% variance"),
    colour = "Group"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.35),
    axis.title = element_text(size = 13, colour = "black"),
    axis.text  = element_text(size = 9, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    legend.key.width  = unit(1.1, "cm"),
    legend.key.height = unit(0.55, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

ggsave("PCA_PC3_PC4_by_group.pdf", p_pc34, width = 7.2, height = 4.8, units = "in", device = cairo_pdf, bg = "white")

####Sample correlation heatmap####

ann_col <- data.frame(
  Group     = factor(sample_meta$Group, levels = c("NTC", "siARID5B")),
  Treatment = factor(sample_meta$Treatment, levels = treatments),
  Timepoint = factor(sample_meta$Timepoint, levels = timepoints),
  row.names = sample_meta$Sample
)

ann_colors <- list(
  Group     = c("NTC" = "#0072B2", "siARID5B" = "#D55E00"),
  Treatment = c("Basal" = "#999999", "E2" = "#CC79A7", "MPA" = "#009E73", "cAMP" = "#F0E442"),
  Timepoint = c("72h" = "#deebf7", "96h" = "#9ecae1", "120h" = "#4292c6", "192h" = "#084594")
)

sample_cor <- cor(logCPM, method = "pearson", use = "pairwise.complete.obs")
ann_col <- ann_col[colnames(sample_cor), , drop = FALSE]

sample_hc <- hclust(as.dist(1 - sample_cor), method = "average")

cor_breaks  <- seq(0.90, 1.00, length.out = 101)
cor_colours <- colorRampPalette(c("#2166AC", "#F7F7F7", "#FDAE61", "#D73027"))(100)

p_heatmap <- pheatmap::pheatmap(
  sample_cor,
  color = cor_colours, breaks = cor_breaks,
  cluster_rows = sample_hc, cluster_cols = sample_hc,
  annotation_row = ann_col, annotation_col = ann_col, annotation_colors = ann_colors,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize = 8, fontsize_row = 6, fontsize_col = 6,
  angle_col = 90, treeheight_row = 45, treeheight_col = 45,
  border_color = NA, silent = TRUE
)

p_heatmap_gg <- ggplotify::as.ggplot(p_heatmap)

ggsave("Supplementary_SampleCorrelationHeatmap.pdf", p_heatmap_gg,
       width = 13, height = 11, units = "in", device = cairo_pdf)


#### DEG bar plot, violin plot, UpSet plots####

PADJ_CUT  <- 0.05
LOGFC_CUT <- 0.5

direction_colours <- c("Down" = "#2166AC", "Up" = "#B2182B")

deg_counts_bar <- all_filtered_res %>%
  filter(!is.na(adj.P.Val), adj.P.Val < PADJ_CUT, abs(logFC) > LOGFC_CUT) %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%
  count(Treatment, Timepoint, Direction, name = "Count") %>%
  complete(Treatment = treatments, Timepoint = timepoints, Direction = c("Up", "Down"), fill = list(Count = 0)) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatments),
    Timepoint = factor(Timepoint, levels = timepoints),
    Direction = factor(Direction, levels = c("Down", "Up")),
    SignedCount = ifelse(Direction == "Down", -Count, Count)
  )

p_bar <- ggplot(deg_counts_bar, aes(x = Timepoint, y = SignedCount, fill = Direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.7, color = "black") +
  geom_text(aes(label = ifelse(Count == 0, "", Count), vjust = ifelse(Direction == "Up", -0.35, 1.35)),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = direction_colours, breaks = c("Up", "Down")) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
  facet_wrap(~ Treatment, nrow = 1) +
  labs(x = "Timepoint", y = "Number of DEGs", fill = "Direction") +
  theme_bw(base_size = 18) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text.x  = element_text(size = 16, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.7),
    strip.text = element_text(size = 19, face = "bold", color = "black", margin = margin(t = 7, b = 7)),
    legend.position = "top",
    legend.title = element_text(size = 17, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(0.7, "cm"),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    plot.margin = margin(8, 15, 10, 10)
  )

ggsave("DEG_counts_barplot.pdf", p_bar, width = 12, height = 5.5, units = "in", device = cairo_pdf, bg = "white")

violin_dat <- all_filtered_res %>%
  filter(!is.na(adj.P.Val), adj.P.Val < PADJ_CUT, abs(logFC) > LOGFC_CUT) %>%
  mutate(Direction = factor(ifelse(logFC > LOGFC_CUT, "Up", "Down"), levels = c("Down", "Up")))

library(scales)

GAP <- 0.5; GAP_FACTOR <- 0.12; MID <- 1.5; OUTER_FACTOR <- 0.25

squish_axis_trans <- function(gap = GAP, gap_factor = GAP_FACTOR, mid = MID, outer_factor = OUTER_FACTOR) {
  y_gap_end <- gap * gap_factor
  y_mid_end <- y_gap_end + (mid - gap)
  trans <- function(x) {
    ax <- abs(x)
    y <- ifelse(ax <= gap, ax * gap_factor,
                ifelse(ax <= mid, y_gap_end + (ax - gap), y_mid_end + (ax - mid) * outer_factor))
    sign(x) * y
  }
  inv <- function(y) {
    ay <- abs(y)
    x <- ifelse(ay <= y_gap_end, ay / gap_factor,
                ifelse(ay <= y_mid_end, gap + (ay - y_gap_end), mid + (ay - y_mid_end) / outer_factor))
    sign(y) * x
  }
  trans_new("squish_axis", trans, inv)
}

p_violin <- ggplot(violin_dat, aes(x = Timepoint, y = logFC, fill = Direction, group = interaction(Timepoint, Direction))) +
  geom_violin(width = 0.88, trim = TRUE, scale = "width", alpha = 0.38, colour = NA, position = position_identity()) +
  geom_boxplot(width = 0.20, outlier.shape = NA, colour = "black", linewidth = 0.75, fatten = 2.5, position = position_identity()) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", colour = "grey55", linewidth = 0.6) +
  scale_fill_manual(values = direction_colours, breaks = c("Up", "Down"), drop = FALSE) +
  scale_y_continuous(breaks = c(-6, -5, -4, -3, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 3, 4, 5), expand = expansion(mult = c(0.02, 0.02))) +
  coord_trans(y = squish_axis_trans(), ylim = c(-6, 5), clip = "off") +
  facet_wrap(~ Treatment, nrow = 1) +
  labs(x = "Timepoint", y = expression(bold(log[2] ~ "fold change")), fill = "Direction") +
  theme_bw(base_size = 18) +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 22, face = "bold", colour = "black", margin = margin(r = 10)),
    axis.text.x  = element_text(size = 18, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 17, colour = "black"),
    strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.8),
    strip.text = element_text(size = 20, face = "bold", colour = "black", margin = margin(t = 7, b = 7)),
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 17),
    legend.key.size = grid::unit(0.75, "cm"),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.7),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 10)
  )

ggsave("ARID5B_KD_logFC_violin.pdf", p_violin, width = 14, height = 9, units = "in", device = cairo_pdf, bg = "white")

#####Centroid distance (uses v$E directly) ####

expr <- v$E
N_TOP_VAR_GENES <- 2000
gene_var <- apply(expr, 1, var)
top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:N_TOP_VAR_GENES]
expr <- expr[top_var_genes, ]

centroid_dist_results <- data.frame()
for (tr in treatments) {
  for (tp in timepoints) {
    meta_sub <- sample_meta %>% filter(Treatment == tr, Timepoint == tp)
    expr_sub <- expr[, meta_sub$Sample]
    ntc_cols <- meta_sub$Sample[meta_sub$Group == "NTC"]
    kd_cols  <- meta_sub$Sample[meta_sub$Group == "siARID5B"]
    centroid_ntc <- rowMeans(expr_sub[, ntc_cols])
    centroid_kd  <- rowMeans(expr_sub[, kd_cols])
    distance <- sqrt(sum((centroid_ntc - centroid_kd)^2))
    centroid_dist_results <- rbind(centroid_dist_results, data.frame(Treatment = tr, Timepoint = tp, Distance = distance))
  }
}
centroid_dist_results$Timepoint <- factor(
  centroid_dist_results$Timepoint,
  levels = c("72h", "96h", "120h", "192h")
)
centroid_dist_results$Treatment <- factor(
  centroid_dist_results$Treatment,
  levels = c("Basal", "E2", "MPA", "cAMP")
)

p_centroid_dist <- ggplot(centroid_dist_results, aes(Timepoint, Distance, group = Treatment, color = Treatment)) +
  geom_line(size = 1) + geom_point(size = 3) +
  scale_colour_manual(values = treatment_colours) +
  theme_classic(base_size = 16) +
  labs(y = "siARID5B vs NTC centroid distance", x = NULL)

ggsave("centroid_distance.pdf", p_centroid_dist, width = 8, height = 5, units = "in", device = cairo_pdf, bg = "white")

#### UpSet plots (uses all_filtered_res)####

get_deg_list <- function(sheet_name, direction = c("all", "up", "down"), adj_p_cutoff = 0.05, logfc_cutoff = 0.5) {
  direction <- match.arg(direction)
  res <- all_filtered_res %>% filter(Contrast == sheet_name)
  keep <- switch(direction,
                 all  = res$adj.P.Val < adj_p_cutoff & abs(res$logFC) > logfc_cutoff,
                 up   = res$adj.P.Val < adj_p_cutoff & res$logFC > logfc_cutoff,
                 down = res$adj.P.Val < adj_p_cutoff & res$logFC < -logfc_cutoff
  )
  res$ENSEMBL[keep & !is.na(keep)]
}

deg_lists_all  <- set_names(map(contrast_sheets, get_deg_list, direction = "all"), contrast_sheets)
deg_lists_up   <- set_names(map(contrast_sheets, get_deg_list, direction = "up"), contrast_sheets)
deg_lists_down <- set_names(map(contrast_sheets, get_deg_list, direction = "down"), contrast_sheets)

plot_upset_treatment <- function(deg_lists, treatment, direction_label, base_font_size = 18) {
  timepoints_ordered <- c("72h", "96h", "120h", "192h")
  treat_contrasts <- grep(paste0("_", treatment, "_"), names(deg_lists), value = TRUE)
  deg_lists_treat <- deg_lists[treat_contrasts]
  names(deg_lists_treat) <- sub(paste0("^siARID5B_vs_NTC_", treatment, "_"), "", names(deg_lists_treat))
  available_timepoints <- intersect(timepoints_ordered, names(deg_lists_treat))
  deg_lists_treat <- deg_lists_treat[available_timepoints]
  
  all_genes <- unique(unlist(deg_lists_treat, use.names = FALSE))
  mat <- map_dfc(deg_lists_treat, ~ all_genes %in% .x)
  names(mat) <- available_timepoints
  rownames(mat) <- all_genes
  mat <- mat %>% mutate(across(all_of(available_timepoints), as.logical))
  
  bar_colour <- switch(direction_label, "Up" = "#B2182B", "Down" = "#2166AC", "All" = "#555555")
  
  p <- upset(
    mat, intersect = rev(available_timepoints),
    base_annotations = list(
      "Intersection size" = intersection_size(counts = TRUE, text = list(size = 5), fill = bar_colour) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme(axis.title.y = element_text(size = base_font_size, face = "bold", margin = margin(r = 8)),
              axis.text.y  = element_text(size = base_font_size - 2, colour = "black"))
    ),
    set_sizes = upset_set_size(geom = geom_bar(fill = "#444444", width = 0.65)) +
      theme(axis.title.x = element_text(size = base_font_size, face = "bold", margin = margin(t = 8)),
            axis.text.x  = element_text(size = base_font_size - 2, colour = "black")),
    sort_intersections_by = "cardinality", sort_sets = FALSE, keep_empty_groups = FALSE,
    width_ratio = 0.25, name = ""
  ) +
    labs(title = paste0(treatment, " - ", direction_label, " DEGs")) +
    theme_bw(base_size = base_font_size) +
    theme(
      plot.title = element_text(size = base_font_size + 4, face = "bold", hjust = 0, margin = margin(b = 10)),
      axis.text = element_text(size = base_font_size - 1, colour = "black"),
      axis.title = element_text(size = base_font_size, face = "bold"),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(t = 12, r = 18, b = 12, l = 12)
    )
  p
}

deg_list_groups <- list(All = deg_lists_all, Up = deg_lists_up, Down = deg_lists_down)
output_directory <- "UpSet_plots"
if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

for (treat in treatments) {
  for (direction in names(deg_list_groups)) {
    p_up <- plot_upset_treatment(deg_list_groups[[direction]], treat, direction, base_font_size = 18)
    file_stem <- paste0("upset_siARID5B_vs_NTC_", treat, "_", direction)
    ggsave(file.path(output_directory, paste0(file_stem, ".pdf")), p_up,
           width = 11, height = 7, units = "in", device = cairo_pdf, bg = "white")
  }
}


####ARID5B knockdown efficiency trajectory (uses v$E / sample_meta directly)####

arid5b_ensembl <- "ENSG00000150347"
arid5b_row <- which(rownames(v$E) == arid5b_ensembl)

arid5b_expr <- data.frame(
  Sample = colnames(v$E),
  log2CPM = v$E[arid5b_row, ],
  Group = sample_meta$Group,
  Treatment = factor(sample_meta$Treatment, levels = treatments),
  Timepoint = factor(sample_meta$Timepoint, levels = timepoints),
  Replicate = sample_meta$Replicate
)

p_arid5b_expr <- ggplot(arid5b_expr, aes(x = Timepoint, y = log2CPM, color = Group, group = Group)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  geom_jitter(width = 0.05, alpha = 0.4, size = 1.5) +
  facet_wrap(~ Treatment, nrow = 1) +
  scale_color_manual(values = c("NTC" = "#0072B2", "siARID5B" = "#D55E00")) +
  labs(y = "ARID5B log2-CPM", x = NULL, color = "Group") +
  theme_bw(base_size = 18) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text.x  = element_text(size = 16, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.7),
    strip.text = element_text(size = 19, face = "bold", color = "black", margin = margin(t = 7, b = 7)),
    legend.position = "top",
    legend.title = element_text(size = 17, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.key.size = grid::unit(0.7, "cm"),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    plot.margin = margin(8, 15, 10, 10)
  )

ggsave("ARID5B_knockdown_trajectory.pdf", p_arid5b_expr, width = 12, height = 5, device = cairo_pdf, bg = "white")


####GSEA Hallmark (builds ranked gene lists FROM all_filtered_res)####

padj_gsea      <- 0.05
nes_cap        <- 3
nes_threshold  <- 1.5
min_timepoints <- 2   # min # of timepoints a pathway must pass |NES| & FDR threshold (direction not considered)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

gsea_results <- list()
set.seed(123)

for (sheet_name in contrast_sheets) {
  
  res <- all_filtered_res %>%
    filter(Contrast == sheet_name, is.finite(t)) %>%
    arrange(desc(abs(t))) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  
  gene_list <- res$t
  names(gene_list) <- res$SYMBOL
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_results[[sheet_name]] <- GSEA(
    geneList = gene_list, TERM2GENE = hallmark,
    pvalueCutoff = 1, pAdjustMethod = "BH",
    minGSSize = 10, maxGSSize = 500,
    verbose = FALSE, seed = TRUE
  )
}

gsea_summary <- map_dfr(names(gsea_results), function(x) {
  df <- as.data.frame(gsea_results[[x]])
  df %>%
    mutate(
      Contrast  = x,
      Treatment = sub("siARID5B_vs_NTC_([^_]+)_.*", "\\1", x),
      Timepoint = sub("siARID5B_vs_NTC_[^_]+_(.*)", "\\1", x),
      pathway   = gsub("HALLMARK_", "", ID),
      pathway   = gsub("_", " ", pathway),
      pathway   = stringr::str_to_title(pathway),
      n_leading_edge = ifelse(is.na(core_enrichment), NA_integer_, stringr::str_count(core_enrichment, "/") + 1)
    ) %>%
    dplyr::select(Contrast, Treatment, Timepoint, pathway, setSize, enrichmentScore, NES,
                  pvalue, p.adjust, qvalue, rank, n_leading_edge, core_enrichment)
})

write.csv(gsea_summary, "All_Contrasts_Hallmark_GSEA_Summary_filtered.csv", row.names = FALSE)
saveRDS(gsea_summary, "All_Contrasts_Hallmark_GSEA_Summary_filtered.rds")

plot_dat <- gsea_summary %>%
  mutate(Treatment = factor(Treatment, levels = treatments), Timepoint = factor(Timepoint, levels = timepoints)) %>%
  complete(pathway, Treatment, Timepoint, fill = list(NES = 0, p.adjust = 1, pvalue = 1)) %>%
  mutate(
    NES_capped     = pmax(pmin(NES, nes_cap), -nes_cap),
    neg_log10_padj = pmin(-log10(pmax(p.adjust, 1e-20)), 20),
    sig            = !is.na(p.adjust) & p.adjust < padj_gsea
  )

pathway_order_all <- plot_dat %>%
  group_by(pathway) %>% summarise(score = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  mutate(score = ifelse(is.finite(score), score, 0)) %>% arrange(score) %>% pull(pathway)

plot_dat_all <- plot_dat %>% mutate(pathway = factor(pathway, levels = pathway_order_all))

hallmark_panel_all <- ggplot(plot_dat_all, aes(x = Timepoint, y = pathway)) +
  geom_point(aes(fill = NES_capped, size = neg_log10_padj), shape = 21, colour = "grey85", stroke = 0.4, na.rm = TRUE) +
  geom_point(data = ~ filter(.x, sig), aes(fill = NES_capped, size = neg_log10_padj),
             shape = 21, colour = "black", stroke = 0.9, na.rm = TRUE) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-nes_cap, nes_cap), breaks = c(-nes_cap, 0, nes_cap),
                       labels = c(paste0("\u2264-", nes_cap), "0", paste0("\u2265", nes_cap)), name = "NES") +
  scale_size_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20), range = c(1.5, 7), name = "-log10(FDR)") +
  labs(
    title = "Hallmark pathway enrichment - siARID5B vs NTC (filtered)",
    subtitle = paste0("All 50 MSigDB hallmarks; black outline = FDR < ", padj_gsea,
                      ". Ranked by maximum |NES| across all contrasts."),
    x = "Timepoint", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 9),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey95"),
    legend.position = "right"
  )

ggsave("Hallmark_panel_ALL_pathways_no_filter_FILTERED_DEGs.pdf", hallmark_panel_all,
       width = 4 + 2.4 * n_distinct(plot_dat_all$Treatment), height = 14, device = cairo_pdf, limitsize = FALSE)

#build_pathway_panels
# Qualifying pathway = passes |NES| >= nes_threshold & FDR < padj_gsea at
# >= min_timepoints timepoints (direction not considered).

build_pathway_panels <- function(plot_dat, top_n_per_treatment = NULL, file_prefix, label_text) {
  
  qualifying_with_score <- plot_dat %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    mutate(pass = abs(NES) >= nes_threshold & p.adjust < padj_gsea) %>%
    filter(pass) %>%
    group_by(Treatment, pathway) %>%
    summarise(
      n_pass      = n_distinct(Timepoint),
      max_abs_NES = max(abs(NES)),
      direction   = ifelse(NES[which.max(abs(NES))] >= 0, "Up", "Down"),  # direction at the ranking timepoint, for labeling only
      .groups = "drop"
    ) %>%
    filter(n_pass >= min_timepoints)
  
  selected <- if (!is.null(top_n_per_treatment)) {
    qualifying_with_score %>% group_by(Treatment) %>%
      slice_max(order_by = max_abs_NES, n = top_n_per_treatment) %>% ungroup()
  } else qualifying_with_score
  
  write.csv(selected %>% arrange(Treatment, desc(max_abs_NES)),
            paste0(file_prefix, "_Pathways_per_Treatment.csv"), row.names = FALSE)
  
  pathway_levels <- selected %>% arrange(Treatment, max_abs_NES) %>%
    mutate(pathway_plot = paste(as.character(Treatment), pathway, sep = "___")) %>%
    pull(pathway_plot) %>% unique()
  
  plot_dat_selected <- plot_dat %>%
    inner_join(selected %>% select(Treatment, pathway, direction), by = c("Treatment", "pathway")) %>%
    mutate(pathway_plot = factor(paste(as.character(Treatment), pathway, sep = "___"), levels = pathway_levels))
  
  p_per_treatment <- ggplot(plot_dat_selected, aes(x = Timepoint, y = pathway_plot)) +
    geom_point(aes(fill = NES_capped, size = neg_log10_padj), shape = 21, colour = "grey85", stroke = 0.4, na.rm = TRUE) +
    geom_point(data = ~ filter(.x, sig), aes(fill = NES_capped, size = neg_log10_padj),
               shape = 21, colour = "black", stroke = 0.9, na.rm = TRUE) +
    facet_wrap(~Treatment, nrow = 1, scales = "free_y") +
    scale_y_discrete(labels = function(x) str_wrap(sub("^.*___", "", x), width = 20)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-nes_cap, nes_cap), breaks = c(-nes_cap, 0, nes_cap),
                         labels = c(paste0("\u2264-", nes_cap), "0", paste0("\u2265", nes_cap)), name = "NES") +
    scale_size_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20), range = c(1.5, 7), name = "-log10(FDR)") +
    labs(
      title = paste0(label_text, " Hallmark pathways per treatment - siARID5B vs NTC"),
      subtitle = paste0("|NES| \u2265 ", nes_threshold, ", FDR < ", padj_gsea, " in \u2265", min_timepoints, " timepoints",
                        if (!is.null(top_n_per_treatment)) paste0("; top ", top_n_per_treatment, " ranked by max |NES|") else "; all qualifying pathways"),
      x = "Timepoint", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 8, lineheight = 0.85), axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey95"),
      legend.position = "right"
    )
  
  n_max <- plot_dat_selected %>% group_by(Treatment) %>% summarise(n = n_distinct(pathway), .groups = "drop") %>% pull(n) %>% max()
  
  ggsave(paste0(file_prefix, "_per_treatment.pdf"), p_per_treatment,
         width = 4 + 2.4 * length(treatments), height = max(4, 0.3 * n_max + 2), device = cairo_pdf, limitsize = FALSE)
  
  union_pathways <- selected %>% distinct(pathway) %>% pull(pathway) %>% as.character()
  plot_dat_combined <- plot_dat %>% filter(as.character(pathway) %in% union_pathways) %>% mutate(pathway = as.character(pathway))
  pathway_order_combined <- plot_dat_combined %>% group_by(pathway) %>%
    summarise(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>% arrange(max_abs_NES) %>% pull(pathway)
  plot_dat_combined <- plot_dat_combined %>% mutate(pathway = factor(pathway, levels = pathway_order_combined))
  
  p_combined <- ggplot(plot_dat_combined, aes(x = Timepoint, y = pathway)) +
    geom_point(aes(fill = NES_capped, size = neg_log10_padj), shape = 21, colour = "grey85", stroke = 0.4, na.rm = TRUE) +
    geom_point(data = ~ filter(.x, sig), aes(fill = NES_capped, size = neg_log10_padj),
               shape = 21, colour = "black", stroke = 0.9, na.rm = TRUE) +
    facet_grid(~ Treatment) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-nes_cap, nes_cap), breaks = c(-nes_cap, 0, nes_cap),
                         labels = c(paste0("\u2264-", nes_cap), "0", paste0("\u2265", nes_cap)), name = "NES") +
    scale_size_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20), range = c(1.5, 7), name = "-log10(FDR)") +
    labs(
      title = "Hallmark pathways - siARID5B vs NTC (union across treatments)",
      subtitle = paste0("|NES| \u2265 ", nes_threshold, ", FDR < ", padj_gsea, " in \u2265", min_timepoints, " timepoints; ", label_text, " pathways per treatment, union across treatments"),
      x = "Timepoint", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 8, lineheight = 0.85), axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey95"),
      panel.spacing = grid::unit(0.6, "lines"), legend.position = "right"
    )
  
  ggsave(paste0(file_prefix, "_COMBINED_union.pdf"), p_combined,
         width = 6 + 1.8 * length(treatments), height = max(6, 0.3 * length(union_pathways) + 2),
         device = cairo_pdf, limitsize = FALSE)
  
  invisible(list(selected = selected, plot_dat_selected = plot_dat_selected))
}

result_top5  <- build_pathway_panels(plot_dat, 5, "Hallmark_panel_TOP5", "Top 5")
result_top10 <- build_pathway_panels(plot_dat, 10, "Hallmark_panel_TOP10", "Top 10")
result_all_qualifying <- build_pathway_panels(plot_dat, NULL, "Hallmark_panel_ALL_qualifying", "All qualifying")


####Leading-edge genes loop + UpSet overlap####

output_dir_base <- "Leading edge genes"
if (!dir.exists(output_dir_base)) dir.create(output_dir_base, recursive = TRUE)

min_timepoints_seen <- 2   # min # of leading-edge timepoints a gene must hit (direction not required to match)
top_n <- 30
pathway_fdr_cutoff <- 0.05
pathway_nes_cutoff <- 1.5
gene_fdr_cutoff <- 0.05
lfc_cap <- 3

reorder_within <- function(x, by, within, sep = "___") {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) sub(paste0(sep, ".*$"), "", x), ...)
}

gsea_all <- gsea_summary

gene_data <- all_filtered_res %>%
  filter(is.finite(t), is.finite(logFC)) %>%
  arrange(desc(abs(t))) %>%
  distinct(SYMBOL, Treatment, Timepoint, .keep_all = TRUE) %>%
  dplyr::select(SYMBOL, logFC, t, adj.P.Val, Treatment, Timepoint)

pathway_specs <- tribble(
  ~target_pathway,                                  ~plot_title,                          ~file_prefix,
  "HALLMARK_E2F_TARGETS",                            "E2F TARGETS",                        "E2F_Targets",
  "HALLMARK_G2M_CHECKPOINT",                          "G2M CHECKPOINT",                     "G2M_Checkpoint",
  "HALLMARK_MITOTIC_SPINDLE",                         "MITOTIC SPINDLE",                    "Mitotic_Spindle",
  "HALLMARK_MYC_TARGETS_V1",                          "MYC TARGETS V1",                     "Myc_Targets_V1",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",       "EPITHELIAL MESENCHYMAL TRANSITION",  "Epithelial_Mesenchymal_Transition",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",               "INTERFERON ALPHA RESPONSE",          "Interferon_Alpha_Response",
  "HALLMARK_TGF_BETA_SIGNALING",                      "TGF BETA SIGNALING",                 "TGF_Beta_Signaling",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",                 "CHOLESTEROL HOMEOSTASIS",            "Cholesterol_Homeostasis",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",                 "ESTROGEN RESPONSE EARLY",            "Estrogen_Response_Early",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",                 "TNFA SIGNALING VIA NFKB",            "TNFa_Signaling_via_NFkB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",                 "IL6 JAK STAT3 SIGNALING",            "IL6_JAK_STAT3_Signaling",
  "HALLMARK_ANGIOGENESIS",                            "ANGIOGENESIS",                       "Angiogenesis"
)

upset_pathway_prefixes <- c("E2F_Targets", "G2M_Checkpoint", "Mitotic_Spindle", "Myc_Targets_V1")
upset_display_names <- c(
  "E2F_Targets" = "E2F", "G2M_Checkpoint" = "G2M",
  "Mitotic_Spindle" = "Mitotic Spindle", "Myc_Targets_V1" = "MYC V1"
)

## --- run_pathway_leading_edge -----------------------------------------------
####GSEA Hallmark (builds ranked gene lists FROM all_filtered_res)####

padj_gsea      <- 0.05
nes_cap        <- 3
nes_threshold  <- 1.5
min_timepoints <- 2   # min # of timepoints a pathway must pass |NES| & FDR threshold (direction not considered)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

gsea_results <- list()
set.seed(123)

for (sheet_name in contrast_sheets) {
  
  res <- all_filtered_res %>%
    filter(Contrast == sheet_name, is.finite(t)) %>%
    arrange(desc(abs(t))) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  
  gene_list <- res$t
  names(gene_list) <- res$SYMBOL
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_results[[sheet_name]] <- GSEA(
    geneList = gene_list, TERM2GENE = hallmark,
    pvalueCutoff = 1, pAdjustMethod = "BH",
    minGSSize = 10, maxGSSize = 500,
    verbose = FALSE, seed = TRUE
  )
}

gsea_summary <- map_dfr(names(gsea_results), function(x) {
  df <- as.data.frame(gsea_results[[x]])
  df %>%
    mutate(
      Contrast  = x,
      Treatment = sub("siARID5B_vs_NTC_([^_]+)_.*", "\\1", x),
      Timepoint = sub("siARID5B_vs_NTC_[^_]+_(.*)", "\\1", x),
      pathway   = gsub("HALLMARK_", "", ID),
      pathway   = gsub("_", " ", pathway),
      pathway   = stringr::str_to_title(pathway),
      n_leading_edge = ifelse(is.na(core_enrichment), NA_integer_, stringr::str_count(core_enrichment, "/") + 1)
    ) %>%
    dplyr::select(Contrast, Treatment, Timepoint, pathway, setSize, enrichmentScore, NES,
                  pvalue, p.adjust, qvalue, rank, n_leading_edge, core_enrichment)
})

write.csv(gsea_summary, "All_Contrasts_Hallmark_GSEA_Summary_filtered.csv", row.names = FALSE)
saveRDS(gsea_summary, "All_Contrasts_Hallmark_GSEA_Summary_filtered.rds")

plot_dat <- gsea_summary %>%
  mutate(Treatment = factor(Treatment, levels = treatments), Timepoint = factor(Timepoint, levels = timepoints)) %>%
  complete(pathway, Treatment, Timepoint, fill = list(NES = 0, p.adjust = 1, pvalue = 1)) %>%
  mutate(
    NES_capped     = pmax(pmin(NES, nes_cap), -nes_cap),
    neg_log10_padj = pmin(-log10(pmax(p.adjust, 1e-20)), 20),
    sig            = !is.na(p.adjust) & p.adjust < padj_gsea
  )

pathway_order_all <- plot_dat %>%
  group_by(pathway) %>% summarise(score = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  mutate(score = ifelse(is.finite(score), score, 0)) %>% arrange(score) %>% pull(pathway)

plot_dat_all <- plot_dat %>% mutate(pathway = factor(pathway, levels = pathway_order_all))

hallmark_panel_all <- ggplot(plot_dat_all, aes(x = Timepoint, y = pathway)) +
  geom_point(aes(fill = NES_capped, size = neg_log10_padj), shape = 21, colour = "grey85", stroke = 0.4, na.rm = TRUE) +
  geom_point(data = ~ filter(.x, sig), aes(fill = NES_capped, size = neg_log10_padj),
             shape = 21, colour = "black", stroke = 0.9, na.rm = TRUE) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-nes_cap, nes_cap), breaks = c(-nes_cap, 0, nes_cap),
                       labels = c(paste0("\u2264-", nes_cap), "0", paste0("\u2265", nes_cap)), name = "NES") +
  scale_size_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20), range = c(1.5, 7), name = "-log10(FDR)") +
  labs(
    title = "Hallmark pathway enrichment - siARID5B vs NTC (filtered)",
    subtitle = paste0("All 50 MSigDB hallmarks; black outline = FDR < ", padj_gsea,
                      ". Ranked by maximum |NES| across all contrasts."),
    x = "Timepoint", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 9),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey95"),
    legend.position = "right"
  )

ggsave("Hallmark_panel_ALL_pathways_no_filter_FILTERED_DEGs.pdf", hallmark_panel_all,
       width = 4 + 2.4 * n_distinct(plot_dat_all$Treatment), height = 14, device = cairo_pdf, limitsize = FALSE)

## --- build_pathway_panels ---------------------------------------------------
## Qualifying pathway = passes |NES| >= nes_threshold & FDR < padj_gsea at
## >= min_timepoints timepoints (direction not considered).
## Top-N selection is a SINGLE ranked list per treatment by max |NES| among
## qualifying pathways.

build_pathway_panels <- function(plot_dat, top_n_per_treatment = NULL, file_prefix, label_text) {
  
  qualifying_with_score <- plot_dat %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    mutate(pass = abs(NES) >= nes_threshold & p.adjust < padj_gsea) %>%
    filter(pass) %>%
    group_by(Treatment, pathway) %>%
    summarise(
      n_pass      = n_distinct(Timepoint),
      max_abs_NES = max(abs(NES)),
      direction   = ifelse(NES[which.max(abs(NES))] >= 0, "Up", "Down"),  # direction at the ranking timepoint, for labeling only
      .groups = "drop"
    ) %>%
    filter(n_pass >= min_timepoints)
  
  selected <- if (!is.null(top_n_per_treatment)) {
    qualifying_with_score %>% group_by(Treatment) %>%
      slice_max(order_by = max_abs_NES, n = top_n_per_treatment) %>% ungroup()
  } else qualifying_with_score
  
  write.csv(selected %>% arrange(Treatment, desc(max_abs_NES)),
            paste0(file_prefix, "_Pathways_per_Treatment.csv"), row.names = FALSE)
  
  pathway_levels <- selected %>% arrange(Treatment, max_abs_NES) %>%
    mutate(pathway_plot = paste(as.character(Treatment), pathway, sep = "___")) %>%
    pull(pathway_plot) %>% unique()
  
  plot_dat_selected <- plot_dat %>%
    inner_join(selected %>% select(Treatment, pathway, direction), by = c("Treatment", "pathway")) %>%
    mutate(pathway_plot = factor(paste(as.character(Treatment), pathway, sep = "___"), levels = pathway_levels))
  
  p_per_treatment <- ggplot(plot_dat_selected, aes(x = Timepoint, y = pathway_plot)) +
    geom_point(aes(fill = NES_capped, size = neg_log10_padj), shape = 21, colour = "grey85", stroke = 0.4, na.rm = TRUE) +
    geom_point(data = ~ filter(.x, sig), aes(fill = NES_capped, size = neg_log10_padj),
               shape = 21, colour = "black", stroke = 0.9, na.rm = TRUE) +
    facet_wrap(~Treatment, nrow = 1, scales = "free_y") +
    scale_y_discrete(labels = function(x) str_wrap(sub("^.*___", "", x), width = 20)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-nes_cap, nes_cap), breaks = c(-nes_cap, 0, nes_cap),
                         labels = c(paste0("\u2264-", nes_cap), "0", paste0("\u2265", nes_cap)), name = "NES") +
    scale_size_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20), range = c(1.5, 7), name = "-log10(FDR)") +
    labs(
      title = paste0(label_text, " Hallmark pathways per treatment - siARID5B vs NTC"),
      subtitle = paste0("|NES| \u2265 ", nes_threshold, ", FDR < ", padj_gsea, " in \u2265", min_timepoints, " timepoints",
                        if (!is.null(top_n_per_treatment)) paste0("; top ", top_n_per_treatment, " ranked by max |NES|") else "; all qualifying pathways"),
      x = "Timepoint", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 8, lineheight = 0.85), axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey95"),
      legend.position = "right"
    )
  
  n_max <- plot_dat_selected %>% group_by(Treatment) %>% summarise(n = n_distinct(pathway), .groups = "drop") %>% pull(n) %>% max()
  
  ggsave(paste0(file_prefix, "_per_treatment.pdf"), p_per_treatment,
         width = 4 + 2.4 * length(treatments), height = max(4, 0.3 * n_max + 2), device = cairo_pdf, limitsize = FALSE)
  
  union_pathways <- selected %>% distinct(pathway) %>% pull(pathway) %>% as.character()
  plot_dat_combined <- plot_dat %>% filter(as.character(pathway) %in% union_pathways) %>% mutate(pathway = as.character(pathway))
  pathway_order_combined <- plot_dat_combined %>% group_by(pathway) %>%
    summarise(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>% arrange(max_abs_NES) %>% pull(pathway)
  plot_dat_combined <- plot_dat_combined %>% mutate(pathway = factor(pathway, levels = pathway_order_combined))
  
  p_combined <- ggplot(plot_dat_combined, aes(x = Timepoint, y = pathway)) +
    geom_point(aes(fill = NES_capped, size = neg_log10_padj), shape = 21, colour = "grey85", stroke = 0.4, na.rm = TRUE) +
    geom_point(data = ~ filter(.x, sig), aes(fill = NES_capped, size = neg_log10_padj),
               shape = 21, colour = "black", stroke = 0.9, na.rm = TRUE) +
    facet_grid(~ Treatment) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-nes_cap, nes_cap), breaks = c(-nes_cap, 0, nes_cap),
                         labels = c(paste0("\u2264-", nes_cap), "0", paste0("\u2265", nes_cap)), name = "NES") +
    scale_size_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20), range = c(1.5, 7), name = "-log10(FDR)") +
    labs(
      title = "Hallmark pathways - siARID5B vs NTC (union across treatments)",
      subtitle = paste0("|NES| \u2265 ", nes_threshold, ", FDR < ", padj_gsea, " in \u2265", min_timepoints, " timepoints; ", label_text, " pathways per treatment, union across treatments"),
      x = "Timepoint", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 8, lineheight = 0.85), axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey95"),
      panel.spacing = grid::unit(0.6, "lines"), legend.position = "right"
    )
  
  ggsave(paste0(file_prefix, "_COMBINED_union.pdf"), p_combined,
         width = 6 + 1.8 * length(treatments), height = max(6, 0.3 * length(union_pathways) + 2),
         device = cairo_pdf, limitsize = FALSE)
  
  invisible(list(selected = selected, plot_dat_selected = plot_dat_selected))
}

result_top5  <- build_pathway_panels(plot_dat, 5, "Hallmark_panel_TOP5", "Top 5")
result_top10 <- build_pathway_panels(plot_dat, 10, "Hallmark_panel_TOP10", "Top 10")
result_all_qualifying <- build_pathway_panels(plot_dat, NULL, "Hallmark_panel_ALL_qualifying", "All qualifying")


####Leading-edge genes loop + UpSet overlap####


output_dir_base <- "Leading edge genes"
if (!dir.exists(output_dir_base)) dir.create(output_dir_base, recursive = TRUE)

min_timepoints_seen <- 2   # min # of leading-edge timepoints a gene must pass threshold at (direction not required to match)
top_n <- 30
pathway_fdr_cutoff <- 0.05
pathway_nes_cutoff <- 1.5
gene_fdr_cutoff <- 0.05
lfc_cap <- 3

reorder_within <- function(x, by, within, sep = "___") {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) sub(paste0(sep, ".*$"), "", x), ...)
}

gsea_all <- gsea_summary

gene_data <- all_filtered_res %>%
  filter(is.finite(t), is.finite(logFC)) %>%
  arrange(desc(abs(t))) %>%
  distinct(SYMBOL, Treatment, Timepoint, .keep_all = TRUE) %>%
  dplyr::select(SYMBOL, logFC, t, adj.P.Val, Treatment, Timepoint)

pathway_specs <- tribble(
  ~target_pathway,                                  ~plot_title,                          ~file_prefix,
  "HALLMARK_E2F_TARGETS",                            "E2F TARGETS",                        "E2F_Targets",
  "HALLMARK_G2M_CHECKPOINT",                          "G2M CHECKPOINT",                     "G2M_Checkpoint",
  "HALLMARK_MITOTIC_SPINDLE",                         "MITOTIC SPINDLE",                    "Mitotic_Spindle",
  "HALLMARK_MYC_TARGETS_V1",                          "MYC TARGETS V1",                     "Myc_Targets_V1",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",       "EPITHELIAL MESENCHYMAL TRANSITION",  "Epithelial_Mesenchymal_Transition",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",               "INTERFERON ALPHA RESPONSE",          "Interferon_Alpha_Response",
  "HALLMARK_TGF_BETA_SIGNALING",                      "TGF BETA SIGNALING",                 "TGF_Beta_Signaling",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",                 "CHOLESTEROL HOMEOSTASIS",            "Cholesterol_Homeostasis",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",                 "ESTROGEN RESPONSE EARLY",            "Estrogen_Response_Early",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",                 "TNFA SIGNALING VIA NFKB",            "TNFa_Signaling_via_NFkB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",                 "IL6 JAK STAT3 SIGNALING",            "IL6_JAK_STAT3_Signaling",
  "HALLMARK_ANGIOGENESIS",                            "ANGIOGENESIS",                       "Angiogenesis"
)

upset_pathway_prefixes <- c("E2F_Targets", "G2M_Checkpoint", "Mitotic_Spindle", "Myc_Targets_V1")
upset_display_names <- c(
  "E2F_Targets" = "E2F", "G2M_Checkpoint" = "G2M",
  "Mitotic_Spindle" = "Mitotic Spindle", "Myc_Targets_V1" = "MYC V1"
)

#run_pathway_leading_edge 
#filter: adj.P.Val < gene_fdr_cutoff in >= min_timepoints_
#Ranking is by the largest |logFC| among PASSING (significant) timepoints


run_pathway_leading_edge <- function(target_pathway, plot_title, file_prefix) {
  
  output_dir <- file.path(output_dir_base, file_prefix)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  target_pathway_label <- target_pathway %>% str_remove("^HALLMARK_") %>% str_replace_all("_", " ") %>% str_to_title()
  
  target_gsea <- gsea_all %>%
    filter(pathway == target_pathway_label, !is.na(NES), !is.na(p.adjust),
           abs(NES) >= pathway_nes_cutoff, p.adjust < pathway_fdr_cutoff,
           !is.na(core_enrichment), core_enrichment != "")
  
  if (nrow(target_gsea) == 0) return(invisible(list(gene_ranking = NULL, union_genes = character(0))))
  
  leading_edge_hits <- target_gsea %>%
    select(Treatment, Timepoint, core_enrichment) %>%
    separate_rows(core_enrichment, sep = "/") %>%
    transmute(Treatment, Timepoint, SYMBOL = str_trim(core_enrichment)) %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct()
  
  gene_passing_dat <- leading_edge_hits %>%
    inner_join(gene_data, by = c("Treatment", "Timepoint", "SYMBOL")) %>%
    filter(!is.na(t), !is.na(adj.P.Val), adj.P.Val < gene_fdr_cutoff)
  
  gene_ranking_full <- gene_passing_dat %>%
    group_by(Treatment, SYMBOL) %>%
    summarise(
      n_leading_edge_timepoints = n_distinct(Timepoint),
      max_abs_logFC             = max(abs(logFC)),
      mean_abs_t_leading_edge   = mean(abs(t), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_leading_edge_timepoints >= min_timepoints_seen) %>%
    group_by(Treatment) %>%
    arrange(desc(max_abs_logFC), .by_group = TRUE) %>%
    mutate(gene_rank = row_number()) %>% ungroup()
  
  gene_ranking <- gene_ranking_full %>% filter(gene_rank <= top_n)
  
  if (nrow(gene_ranking) == 0) return(invisible(list(gene_ranking = NULL, union_genes = character(0), union_genes_uncapped = character(0))))
  
  write.csv(gene_ranking, file.path(output_dir, paste0(file_prefix, "_Recurrent_LeadingEdge_Top", top_n, ".csv")), row.names = FALSE)
  
  all_treatment_data <- gene_data %>%
    semi_join(gene_ranking, by = c("Treatment", "SYMBOL")) %>%
    left_join(gene_ranking, by = c("Treatment", "SYMBOL"))
  
  plot_dat_le <- all_treatment_data %>%
    mutate(
      Treatment = factor(Treatment, levels = treatments), Timepoint = factor(Timepoint, levels = timepoints),
      logFC_capped = pmax(pmin(logFC, lfc_cap), -lfc_cap),
      sig = !is.na(adj.P.Val) & adj.P.Val < gene_fdr_cutoff,
      SYMBOL_ord = reorder_within(SYMBOL, -gene_rank, Treatment)
    )
  
  p_le <- ggplot(plot_dat_le, aes(x = Timepoint, y = SYMBOL_ord, fill = logFC_capped)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(data = plot_dat_le %>% filter(sig), aes(label = "*"), colour = "black", size = 4.5, fontface = "bold", vjust = 0.75) +
    facet_wrap(~ Treatment, nrow = 1, scales = "free_y") +
    scale_y_reordered() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lfc_cap, lfc_cap), breaks = c(-lfc_cap, 0, lfc_cap),
                         labels = c(paste0("\u2264-", lfc_cap), "0", paste0("\u2265", lfc_cap)), name = "siARID5B vs NTC\nlog2FC") +
    labs(
      title = plot_title,
      subtitle = paste0("Leading-edge genes with FDR < ", gene_fdr_cutoff, " at \u2265", min_timepoints_seen, " timepoints; top ", top_n, " ranked by max |log2FC| among significant timepoints"),
      caption = paste0("Pathway qualifying criterion: |NES| \u2265 ", pathway_nes_cutoff, " and pathway FDR < ", pathway_fdr_cutoff, ". * Gene-level FDR < ", gene_fdr_cutoff, "."),
      x = "Timepoint", y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey30"),
      plot.caption = element_text(hjust = 0, size = 8, colour = "grey35"),
      strip.background = element_rect(fill = "grey85", colour = "grey50"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 6.5), axis.text.x = element_text(size = 8),
      panel.grid = element_blank(), panel.spacing = grid::unit(0.6, "lines"),
      legend.title = element_text(face = "bold", size = 9), legend.text = element_text(size = 8), legend.position = "right"
    )
  
  max_genes <- plot_dat_le %>% group_by(Treatment) %>% summarise(n_genes = n_distinct(SYMBOL), .groups = "drop") %>% pull(n_genes) %>% max()
  
  ggsave(file.path(output_dir, paste0(file_prefix, "_Recurrent_LeadingEdge_Top", top_n, ".pdf")), plot = p_le,
         width = 4 * length(treatments), height = max(6, 0.23 * max_genes + 2), units = "in", device = cairo_pdf, limitsize = FALSE)
  
  union_genes <- gene_ranking %>% distinct(SYMBOL) %>% pull(SYMBOL)
  
  union_gene_order <- leading_edge_hits %>%
    filter(SYMBOL %in% union_genes) %>%
    inner_join(gene_data, by = c("Treatment", "Timepoint", "SYMBOL")) %>%
    group_by(SYMBOL) %>%
    summarise(n_leading_edge_conditions = n_distinct(paste(Treatment, Timepoint)),
              mean_abs_t_leading_edge = mean(abs(t), na.rm = TRUE), .groups = "drop") %>%
    arrange(n_leading_edge_conditions, mean_abs_t_leading_edge)
  
  write.csv(union_gene_order %>% arrange(desc(n_leading_edge_conditions), desc(mean_abs_t_leading_edge)),
            file.path(output_dir, paste0(file_prefix, "_LeadingEdge_Union_GeneRanking.csv")), row.names = FALSE)
  
  plot_dat_union <- gene_data %>%
    filter(SYMBOL %in% union_genes) %>%
    left_join(union_gene_order, by = "SYMBOL") %>%
    mutate(
      Treatment = factor(Treatment, levels = treatments), Timepoint = factor(Timepoint, levels = timepoints),
      SYMBOL = factor(SYMBOL, levels = union_gene_order$SYMBOL),
      logFC_capped = pmax(pmin(logFC, lfc_cap), -lfc_cap),
      sig = !is.na(adj.P.Val) & adj.P.Val < gene_fdr_cutoff
    )
  
  p_union <- ggplot(plot_dat_union, aes(x = Timepoint, y = SYMBOL, fill = logFC_capped)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(data = plot_dat_union %>% filter(sig), aes(label = "*"), colour = "black", size = 4, fontface = "bold", vjust = 0.75) +
    facet_grid(~ Treatment) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-lfc_cap, lfc_cap), breaks = c(-lfc_cap, 0, lfc_cap),
                         labels = c(paste0("\u2264-", lfc_cap), "0", paste0("\u2265", lfc_cap)), name = "siARID5B vs NTC\nlog2FC") +
    labs(
      title = plot_title,
      subtitle = paste0("Union of the top ", top_n, " leading-edge genes selected within each treatment (n = ", length(union_genes), " genes)"),
      caption = paste0("Genes ordered by total leading-edge recurrence, then mean |moderated t|. * Gene-level FDR < ", gene_fdr_cutoff, "."),
      x = "Timepoint", y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey30"),
      plot.caption = element_text(hjust = 0, size = 8, colour = "grey35"),
      strip.background = element_rect(fill = "grey85", colour = "grey50"),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 6.5), axis.text.x = element_text(size = 8),
      panel.grid = element_blank(), panel.spacing = grid::unit(0.5, "lines"),
      legend.title = element_text(face = "bold", size = 9), legend.text = element_text(size = 8), legend.position = "right"
    )
  
  n_union_genes <- length(union_genes)
  ggsave(file.path(output_dir, paste0(file_prefix, "_LeadingEdge_Union_AllTreatments.pdf")), plot = p_union,
         width = 14, height = max(7, 0.22 * n_union_genes + 2), units = "in", device = cairo_pdf, limitsize = FALSE)
  
  union_genes_uncapped <- gene_ranking_full %>% distinct(SYMBOL) %>% pull(SYMBOL)
  
  invisible(list(gene_ranking = gene_ranking, union_genes = union_genes, union_genes_uncapped = union_genes_uncapped))
}

results_le <- pmap(pathway_specs, function(target_pathway, plot_title, file_prefix) {
  run_pathway_leading_edge(target_pathway, plot_title, file_prefix)
})
names(results_le) <- pathway_specs$file_prefix

pathway_gene_sets <- map(upset_pathway_prefixes, ~ results_le[[.x]]$union_genes_uncapped)
names(pathway_gene_sets) <- unname(upset_display_names[upset_pathway_prefixes])

all_genes_upset <- sort(unique(unlist(pathway_gene_sets)))
membership_table <- tibble(SYMBOL = all_genes_upset)
for (nm in names(pathway_gene_sets)) {
  membership_table[[nm]] <- membership_table$SYMBOL %in% pathway_gene_sets[[nm]]
}

membership_table <- membership_table %>%
  rowwise() %>%
  mutate(
    n_pathways = sum(c_across(all_of(names(pathway_gene_sets)))),
    pathways = paste(names(pathway_gene_sets)[c_across(all_of(names(pathway_gene_sets)))], collapse = " & ")
  ) %>%
  ungroup() %>%
  arrange(desc(n_pathways), pathways, SYMBOL)

write.csv(membership_table, file.path(output_dir_base, "Pathway_LeadingEdge_Gene_Membership_filtered.csv"), row.names = FALSE)

p_upset <- ComplexUpset::upset(
  membership_table, intersect = names(pathway_gene_sets), name = "Pathway",
  min_size = 1, width_ratio = 0.15, sort_intersections_by = "cardinality",
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(counts = TRUE, bar_number_threshold = 1, fill = "#4a4a8a") +
      ggplot2::labs(y = "Intersection size\n(# genes)")
  ),
  set_sizes = ComplexUpset::upset_set_size() + ggplot2::labs(y = "Union set size\n(# genes)"),
  themes = ComplexUpset::upset_default_themes(text = ggplot2::element_text(size = 11))
) + ggplot2::labs(title = "Leading-Edge Gene Overlap Across Pathways (no rank cap)")

ggplot2::ggsave(file.path(output_dir_base, "Pathway_LeadingEdge_UpSet_filtered.pdf"), plot = p_upset,
                width = 9, height = 6, units = "in", device = cairo_pdf)



####Cytokine/chemokine survey####

gene_classes <- list(
  "Chemokine ligand" = c(
    "CCL1","CCL2","CCL3","CCL3L1","CCL3L3","CCL4","CCL4L1","CCL4L2",
    "CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15","CCL16",
    "CCL17","CCL18","CCL19","CCL20","CCL21","CCL22","CCL23","CCL24",
    "CCL25","CCL26","CCL27","CCL28",
    "CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","CXCL10",
    "CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17",
    "XCL1","XCL2","CX3CL1"
  ),
  "Chemokine receptor" = c(
    "CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10",
    "CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6",
    "XCR1","CX3CR1",
    "ACKR1","ACKR2","ACKR3","ACKR4","GPR182","PITPNM3","CCRL2"
  ),
  "Interleukin" = c(
    "IL1A","IL1B","IL1RN","IL2","IL3","IL4","IL5","IL6","IL7","IL9",
    "IL10","IL11","IL12A","IL12B","IL13","IL15","IL16","IL17A","IL17B",
    "IL17C","IL17D","IL17F","IL18","IL19","IL20","IL21","IL22","IL23A",
    "IL24","IL25","IL26","IL27","IL31","IL32","IL33","IL34",
    "IL36A","IL36B","IL36G","IL36RN","IL37","IL1F10"
  ),
  "TNF superfamily" = c(
    "TNF","LTA","LTB","TNFSF4","TNFSF8","TNFSF9","TNFSF10","TNFSF11",
    "TNFSF12","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18",
    "EDA","FASLG","CD40LG","CD70"
  ),
  "Colony-stimulating factor" = c("CSF1","CSF2","CSF3"),
  "Interferon" = c(
    "IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10",
    "IFNA13","IFNA14","IFNA16","IFNA17","IFNA21",
    "IFNB1","IFNE","IFNG","IFNK","IFNW1",
    "IFNL1","IFNL2","IFNL3","IFNL4"
  )
)

class_lookup <- imap_dfr(gene_classes, ~ tibble(SYMBOL = .x, Gene_class = .y)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

master_gene_list <- class_lookup$SYMBOL

fdr_cutoff_cytokine <- 0.05
lfc_cap_cytokine    <- 3

chemokine_cytokine_survey <- all_filtered_res %>%
  filter(SYMBOL %in% master_gene_list) %>%
  left_join(class_lookup, by = "SYMBOL") %>%
  mutate(
    neg_log10_FDR = -log10(pmax(adj.P.Val, 1e-10)),
    logFC_capped  = pmax(pmin(logFC, lfc_cap_cytokine), -lfc_cap_cytokine),
    sig           = !is.na(adj.P.Val) & adj.P.Val < fdr_cutoff_cytokine
  )

write.xlsx(chemokine_cytokine_survey, "chemokine_cytokine_survey_filtered.xlsx", overwrite = TRUE)

plot_dat_cytokine <- chemokine_cytokine_survey %>%
  dplyr::select(SYMBOL, Gene_class, Treatment, Timepoint, logFC_capped, adj.P.Val, neg_log10_FDR, sig) %>%
  distinct(SYMBOL, Treatment, Timepoint, .keep_all = TRUE) %>%
  complete(nesting(SYMBOL, Gene_class), Treatment, Timepoint)

p_cytokine <- ggplot(plot_dat_cytokine, aes(x = Timepoint, y = SYMBOL)) +
  geom_point(aes(fill = logFC_capped, size = neg_log10_FDR), shape = 21, colour = "grey85", stroke = 0.4, na.rm = TRUE) +
  geom_point(data = ~ filter(.x, sig), aes(fill = logFC_capped, size = neg_log10_FDR),
             shape = 21, colour = "black", stroke = 0.9, na.rm = TRUE) +
  facet_grid(Gene_class ~ Treatment, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lfc_cap_cytokine, lfc_cap_cytokine), name = "log2FC\n(siARID5B vs NTC)") +
  scale_size_continuous(range = c(1.5, 6), name = "-log10(FDR)") +
  labs(
    title = "Cytokine and chemokine expression - siARID5B vs NTC",
    subtitle = paste0("Black outline = FDR < ", fdr_cutoff_cytokine, "."),
    x = "Timepoint", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 8, colour = "grey40"),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text.x = element_text(face = "bold", size = 11),
    strip.text.y = element_text(angle = 0, size = 9, hjust = 0),
    axis.text.y = element_text(size = 7),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey95"),
    legend.position = "right"
  )

ggsave("cytokine_dotplot_all_classes_filtered.pdf", p_cytokine,
       width = 4 + 2.2 * nlevels(plot_dat_cytokine$Treatment),
       height = max(6, 0.28 * n_distinct(plot_dat_cytokine$SYMBOL) + 0.4 * n_distinct(plot_dat_cytokine$Gene_class) + 2),
       device = cairo_pdf, limitsize = FALSE)


####Leading edge genes dotplot 1####

fdr_cutoff_curated <- 0.05
lfc_cap_curated     <- 3

gene_sets <- list(
  #"TGF-\u03B2 Signaling" = c("TGFBR1", "CDKN1C", "SMAD3", "SERPINE1", "ID2", "ID3"),
  "Progesterone receptor signaling" = c("PGR", "FOSL2", "FKBP5", "FKBP4", "NCOA1"),
  "Decidualization" = c("HOXA10", "PRL", "IGFBP1", "LEFTY2", "FOXO1", "WNT4", "HSD11B1"),
  "Cholesterol biosynthesis" = c("SREBF2", "HMGCR", "HMGCS1", "FDPS", "MVD"),
  "Angiogenesis" = c("OLR1", "SPP1", "VEGFA"),
  "E2F targets / proliferation" = c("EZH2", "TOP2A", "AURKA"),
  "Epithelial-mesenchymal transition" = c("CDH2", "COL11A1", "THBS1"),
  "TNFA-NFKB" = c("TNFAIP2", "BCL3", "STAT5A"),
  "IL-6/JAK-STAT3 signaling" = c("TNFRSF12A", "SOCS1", "CSF1", "SOCS3"),
  "TNF superfamily ligands" = c("TNFSF4", "TNFSF18")
)

gene_set_lookup <- imap_dfr(gene_sets, ~ tibble(GeneSet = .y, SYMBOL = .x))
genes_use <- gene_set_lookup$SYMBOL
symbol_levels_curated <- unique(genes_use)

plot_dat_curated <- all_filtered_res %>%
  filter(SYMBOL %in% genes_use) %>%
  distinct(SYMBOL, Treatment, Timepoint, .keep_all = TRUE) %>%
  dplyr::select(SYMBOL, Treatment, Timepoint, logFC, adj.P.Val) %>%
  inner_join(gene_set_lookup, by = "SYMBOL") %>%
  mutate(
    GeneSet      = factor(GeneSet, levels = names(gene_sets)),
    SYMBOL       = factor(SYMBOL, levels = rev(symbol_levels_curated)),
    Treatment    = factor(Treatment, levels = treatments),
    Timepoint    = factor(Timepoint, levels = timepoints),
    logFC_capped = pmax(pmin(logFC, lfc_cap_curated), -lfc_cap_curated),
    neglog10FDR  = -log10(pmax(adj.P.Val, 1e-10)),
    sig          = !is.na(adj.P.Val) & adj.P.Val < fdr_cutoff_curated
  )

p_curated <- ggplot(plot_dat_curated, aes(x = Timepoint, y = SYMBOL)) +
  geom_point(aes(size = neglog10FDR, fill = logFC_capped), shape = 21, colour = "grey75", stroke = 0.5, na.rm = TRUE) +
  geom_point(data = plot_dat_curated %>% filter(sig), aes(size = neglog10FDR, fill = logFC_capped),
             shape = 21, colour = "black", stroke = 1.3, na.rm = TRUE) +
  facet_grid(GeneSet ~ Treatment, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lfc_cap_curated, lfc_cap_curated), breaks = c(-lfc_cap_curated, 0, lfc_cap_curated),
                       labels = c(paste0("\u2264-", lfc_cap_curated), "0", paste0("\u2265", lfc_cap_curated)), name = "log2FC") +
  scale_size_continuous(range = c(2.5, 8), name = expression(bold(-log[10]("FDR")))) +
  labs(x = NULL, y = NULL, title = "Curated Marker Genes Across Treatment and Time") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.margin = margin(t = 20, r = 10, b = 5, l = 5),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text.x = element_text(face = "bold", size = 14),
    strip.text.y = element_text(face = "bold", size = 11, angle = 0),
    axis.text.y = element_text(size = 11, face = "italic"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey94"), panel.grid.minor = element_blank(),
    panel.spacing = unit(0.4, "lines"),
    legend.title = element_text(face = "bold", size = 13), legend.text = element_text(size = 11),
    legend.box.spacing = unit(0.8, "cm")
  )

n_total_genes_curated <- length(symbol_levels_curated)

ggsave("Curated_GeneSets_Combined_Dotplot.pdf", p_curated,
       width = 12, height = max(6, 0.55 * n_total_genes_curated + 3), units = "in",
       device = cairo_pdf, limitsize = FALSE)


####Leading edge genes dotplot 2####
#take top 4 genes (1 from each treatment/4 from one treatment if other treatments not significant)
fdr_cutoff_curated <- 0.05
lfc_cap_curated     <- 3

gene_sets <- list(
  #"TGF-\u03B2 Signaling" = c("TGFBR1", "CDKN1C", "SMAD3", "SERPINE1", "ID2", "ID3"),
  "Progesterone receptor signaling" = c("PGR", "FOSL2", "FKBP5", "FKBP4", "NCOA1"),
  "Decidualization" = c("HOXA10", "PRL", "IGFBP1", "LEFTY2", "FOXO1", "WNT4", "HSD11B1"),
  "Cholesterol homeostasis" = c("ACAT2", "CD9", "ALCAM"), #1 from basal, 1 from MPA, 1 from cAMP
  "Angiogenesis" = c("POSTN", "SPP1","LUM"), #1 from basal, 1 from E2, 1 from MPA
  "E2F targets / proliferation" = c("E2F8", "DLGAP5", "BIRC5","CDC20"), #1 from basal, cAMP, 1 from E2, 1 from MPA
  "Epithelial-mesenchymal transition" = c("POSTN", "MGP", "SERPINE2","COL11A1"), # 1 from each treatment
  "TNFA-NFKB" = c("BHLHE40", "PLAU"), #1 from MPA, 1 from cAMP
  "IL-6/JAK-STAT3 signaling" = c("CD14", "SOCS3", "ITGB3","SOCS1"), # 3 from cAMP
  "TNF superfamily ligands" = c("TNFSF4", "TNFSF18")
)

gene_set_lookup <- imap_dfr(gene_sets, ~ tibble(GeneSet = .y, SYMBOL = .x))
genes_use <- gene_set_lookup$SYMBOL
symbol_levels_curated <- unique(genes_use)

plot_dat_curated <- all_filtered_res %>%
  filter(SYMBOL %in% genes_use) %>%
  distinct(SYMBOL, Treatment, Timepoint, .keep_all = TRUE) %>%
  dplyr::select(SYMBOL, Treatment, Timepoint, logFC, adj.P.Val) %>%
  inner_join(gene_set_lookup, by = "SYMBOL") %>%
  mutate(
    GeneSet      = factor(GeneSet, levels = names(gene_sets)),
    SYMBOL       = factor(SYMBOL, levels = rev(symbol_levels_curated)),
    Treatment    = factor(Treatment, levels = treatments),
    Timepoint    = factor(Timepoint, levels = timepoints),
    logFC_capped = pmax(pmin(logFC, lfc_cap_curated), -lfc_cap_curated),
    neglog10FDR  = -log10(pmax(adj.P.Val, 1e-10)),
    sig          = !is.na(adj.P.Val) & adj.P.Val < fdr_cutoff_curated
  )

p_curated <- ggplot(plot_dat_curated, aes(x = Timepoint, y = SYMBOL)) +
  geom_point(aes(size = neglog10FDR, fill = logFC_capped), shape = 21, colour = "grey75", stroke = 0.5, na.rm = TRUE) +
  geom_point(data = plot_dat_curated %>% filter(sig), aes(size = neglog10FDR, fill = logFC_capped),
             shape = 21, colour = "black", stroke = 1.3, na.rm = TRUE) +
  facet_grid(GeneSet ~ Treatment, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lfc_cap_curated, lfc_cap_curated), breaks = c(-lfc_cap_curated, 0, lfc_cap_curated),
                       labels = c(paste0("\u2264-", lfc_cap_curated), "0", paste0("\u2265", lfc_cap_curated)), name = "log2FC") +
  scale_size_continuous(range = c(2.5, 8), name = expression(bold(-log[10]("FDR")))) +
  labs(x = NULL, y = NULL, title = "Curated Marker Genes Across Treatment and Time") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.margin = margin(t = 20, r = 10, b = 5, l = 5),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text.x = element_text(face = "bold", size = 14),
    strip.text.y = element_text(face = "bold", size = 11, angle = 0),
    axis.text.y = element_text(size = 11, face = "italic"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey94"), panel.grid.minor = element_blank(),
    panel.spacing = unit(0.4, "lines"),
    legend.title = element_text(face = "bold", size = 13), legend.text = element_text(size = 11),
    legend.box.spacing = unit(0.8, "cm")
  )

n_total_genes_curated <- length(symbol_levels_curated)

ggsave("Curated_GeneSets_Combined_logFC_Dotplot.pdf", p_curated,
       width = 12, height = max(6, 0.55 * n_total_genes_curated + 3), units = "in",
       device = cairo_pdf, limitsize = FALSE)


####PROGENy####

progeny_top_genes  <- 500
progeny_FDR_cutoff <- 0.05

contrast_statistics <- set_names(contrast_sheets) %>%
  map(function(sheet_name) {
    all_filtered_res %>%
      filter(Contrast == sheet_name) %>%
      transmute(SYMBOL = SYMBOL %>% as.character() %>% str_trim() %>% str_to_upper(), t = as.numeric(t)) %>%
      filter(!is.na(t), is.finite(t)) %>%
      arrange(desc(abs(t))) %>%
      distinct(SYMBOL, .keep_all = TRUE) %>%
      tibble::deframe()
  }) %>%
  purrr::compact()

common_genes <- Reduce(intersect, lapply(contrast_statistics, names))

t_statistic_matrix <- do.call(cbind, lapply(contrast_statistics, function(x) x[common_genes]))
rownames(t_statistic_matrix) <- common_genes
colnames(t_statistic_matrix) <- names(contrast_statistics)
storage.mode(t_statistic_matrix) <- "double"

progeny_network <- decoupleR::get_progeny(organism = "human", top = progeny_top_genes)

progeny_results_raw <- decoupleR::run_mlm(
  mat = t_statistic_matrix, network = progeny_network,
  .source = "source", .target = "target", .mor = "weight", minsize = 5
)

progeny_results <- progeny_results_raw %>%
  rename(Pathway = source, Contrast = condition, PROGENy_score = score, PROGENy_pvalue = p_value) %>%
  mutate(
    Contrast  = as.character(Contrast),
    Treatment = str_match(Contrast, "^siARID5B_vs_NTC_([^_]+)_(.+)$")[, 2],
    Timepoint = str_match(Contrast, "^siARID5B_vs_NTC_([^_]+)_(.+)$")[, 3]
  ) %>%
  group_by(Contrast) %>%
  mutate(PROGENy_FDR = stats::p.adjust(PROGENy_pvalue, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Treatment = factor(Treatment, levels = treatments), Timepoint = factor(Timepoint, levels = timepoints),
    Significant = !is.na(PROGENy_FDR) & PROGENy_FDR < progeny_FDR_cutoff,
    neg_log10_FDR = if_else(is.na(PROGENy_FDR), NA_real_, -log10(pmax(PROGENy_FDR, .Machine$double.xmin))),
    Activity_direction = case_when(
      PROGENy_score > 0 ~ "Higher activity in ARID5B KD",
      PROGENy_score < 0 ~ "Lower activity in ARID5B KD",
      TRUE ~ "No difference"
    )
  )

write.csv(progeny_results %>% arrange(Treatment, Timepoint, Pathway),
          "All_Contrasts_PROGENy_MLM_Pathway_Activities_filtered.csv", row.names = FALSE)

activity_color_limit <- progeny_results %>% summarise(limit = quantile(abs(PROGENy_score), probs = 0.95, na.rm = TRUE)) %>% pull(limit) %>% ceiling()
activity_color_limit <- max(3, activity_color_limit)

progeny_results <- progeny_results %>%
  mutate(PROGENy_score_capped = pmax(pmin(PROGENy_score, activity_color_limit), -activity_color_limit))

progeny_pathway_order <- c("Androgen", "Estrogen", "EGFR", "MAPK", "PI3K", "WNT", "TGFb", "JAK-STAT", "NFkB", "TNFa", "Trail", "p53", "Hypoxia", "VEGF")
progeny_results <- progeny_results %>% mutate(Pathway = factor(Pathway, levels = rev(progeny_pathway_order)))

p_progeny_all <- ggplot(progeny_results, aes(x = Timepoint, y = Pathway)) +
  geom_point(aes(fill = PROGENy_score_capped, size = neg_log10_FDR), shape = 21, colour = "grey75", stroke = 0.5, na.rm = TRUE) +
  geom_point(data = progeny_results %>% filter(Significant), aes(fill = PROGENy_score_capped, size = neg_log10_FDR),
             shape = 21, colour = "black", stroke = 1.4, na.rm = TRUE) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_fill_gradient2(name = "PROGENy\nactivity", low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-activity_color_limit, activity_color_limit), breaks = c(-activity_color_limit, 0, activity_color_limit),
                       labels = c(paste0("\u2264-", activity_color_limit), "0", paste0("\u2265", activity_color_limit))) +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(1.5, 8)) +
  labs(title = "PROGENy pathway activity - siARID5B vs NTC (filtered)", x = "Timepoint", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 9),
    panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey95"),
    legend.position = "right"
  )

n_progeny_treatments <- n_distinct(progeny_results$Treatment)
ggsave("PROGENy_All_Pathways_faceted_dotplot_filtered.pdf", p_progeny_all,
       width = 4 + 2.4 * n_progeny_treatments, height = 6, device = cairo_pdf, limitsize = FALSE)



