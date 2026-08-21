setwd("~/Documents/1.Endometriosis/2. Experiment/2. EnS Dfx/4. rt-qPCR/20260804 Cheung 3 batches for paper")
library(tidyverse)
library(readxl)
####main figure####
timepoint_order <- c(
  "iPSC", "D0", "D4", "D8", "D12", "D20", "D28"
)

raw_data <- read_excel(
  "paper_figure.xlsx",
  sheet = "Main figure data"
)

available_timepoints <- timepoint_order[
  timepoint_order %in% colnames(raw_data)
]

# ARID5B on top and CD10 on the bottom
gene_order <- c("ARID5B", "MME")

plot_data <- raw_data %>%
  filter(Gene %in% gene_order) %>%
  pivot_longer(
    cols = all_of(available_timepoints),
    names_to = "Timepoint",
    values_to = "FoldChange"
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = gene_order
    ),
    Timepoint = factor(
      Timepoint,
      levels = available_timepoints
    )
  )

summary_data <- plot_data %>%
  group_by(Gene, Timepoint) %>%
  summarise(
    Mean = mean(FoldChange, na.rm = TRUE),
    SD = sd(FoldChange, na.rm = TRUE),
    .groups = "drop"
  )

p_main <- ggplot() +
  geom_col(
    data = summary_data,
    aes(
      x = Timepoint,
      y = Mean
    ),
    width = 0.70,
    fill = NA,
    color = "black",
    linewidth = 0.65
  ) +
  geom_errorbar(
    data = summary_data,
    aes(
      x = Timepoint,
      ymin = pmax(Mean - SD, 0),
      ymax = Mean + SD
    ),
    width = 0.18,
    linewidth = 0.60,
    color = "black"
  ) +
  geom_point(
    data = plot_data,
    aes(
      x = Timepoint,
      y = FoldChange
    ),
    position = position_jitter(
      width = 0.10,
      height = 0,
      seed = 123
    ),
    shape = 21,
    size = 2.7,
    stroke = 0.65,
    fill = scales::alpha("#6B8E23", 0.5),
    color = "black"
  ) +
  facet_wrap(
    vars(Gene),
    ncol = 1,
    scales = "free_y"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.15)
    )
  ) +
  labs(
    x = NULL,
    y = "Fold change relative to iPSC"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 14,
      color = "black",
      margin = margin(r = 12)
    ),
    
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.line.x = element_line(
      linewidth = 0.65,
      color = "black"
    ),
    
    axis.line.y = element_line(
      linewidth = 0.65,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.50,
      color = "black"
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 16,
      face = "italic",
      color = "black"
    ),
    
    panel.spacing = grid::unit(
      1.3,
      "lines"
    )
  )

p_main

ggsave(
  "main_figure_ARID5B_CD10.pdf",
  plot = p_main,
  width = 4,
  height = 5,
  units = "in",
  device = grDevices::pdf,
  useDingbats = FALSE
)



####supplementary figure####
raw_data <- read_excel(
  "paper_figure.xlsx",
  sheet = "Supp heatmap data"
)

timepoint_order <- c(
  "iPSC", "D0", "D4", "D8", "D12", "D20", "D28"
)



available_timepoints <- timepoint_order[
  timepoint_order %in% colnames(raw_data)
]

gene_order <- unique(raw_data$Gene)

plot_data <- raw_data %>%
  pivot_longer(
    cols = all_of(available_timepoints),
    names_to = "Timepoint",
    values_to = "FoldChange"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Timepoint = factor(
      Timepoint,
      levels = available_timepoints
    )
  )

summary_data <- plot_data %>%
  group_by(Gene, Timepoint) %>%
  summarise(
    Mean = mean(FoldChange, na.rm = TRUE),
    SD = sd(FoldChange, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot() +
  geom_col(
    data = summary_data,
    aes(
      x = Timepoint,
      y = Mean
    ),
    width = 0.70,
    fill = NA,
    color = "black",
    linewidth = 0.65
  ) +
  geom_errorbar(
    data = summary_data,
    aes(
      x = Timepoint,
      ymin = pmax(Mean - SD, 0),
      ymax = Mean + SD
    ),
    width = 0.18,
    linewidth = 0.60,
    color = "black"
  ) +
  geom_point(
    data = plot_data,
    aes(
      x = Timepoint,
      y = FoldChange
    ),
    position = position_jitter(
      width = 0.10,
      height = 0,
      seed = 123
    ),
    shape = 21,
    size = 2.7,
    stroke = 0.65,
    fill = scales::alpha("#6B8E23", 0.5),
    color = "black"
  ) +
  facet_wrap(
    vars(Gene),
    ncol = 4,
    scales = "free_y"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.15)
    )
  ) +
  labs(
    x = NULL,
    y = "Fold change relative to iPSC"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 22,
      color = "black",
      margin = margin(r = 12)
    ),
    
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.line.x = element_line(
      linewidth = 0.65,
      color = "black"
    ),
    
    axis.line.y = element_line(
      linewidth = 0.65,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.50,
      color = "black"
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 14,
      face = "italic",
      color = "black"
    ),
    
    panel.spacing = grid::unit(
      1.1,
      "lines"
    )
  )
p

ggsave(
  "supp_gene_expression_faceted_barplot2.pdf",
  plot = p,
  width = 10,
  height = 10.5,
  units = "in",
  device = grDevices::pdf,
  useDingbats = FALSE
)

