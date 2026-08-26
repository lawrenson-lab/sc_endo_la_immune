#!/usr/bin/env Rscript --vanilla

suppressMessages(library(tidyverse,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/"))
suppressMessages(library(Giotto,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/"))
suppressMessages(library(data.table,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/"))
suppressMessages(library(dendextend,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/"))
suppressMessages(library(igraph,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/"))
suppressMessages(library(factoextra,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/"))


slide<-commandArgs(trailingOnly=TRUE)

start.time <- Sys.time()
counts<-fread(slide)%>%as.data.frame()
counts<-counts%>%mutate(probe=str_replace(probe," ","_"))%>%column_to_rownames("probe")
slide<-str_remove_all(slide,"docs/|asinh.csv")
print(slide)
probes<-rownames(counts)
metadata<-fread("docs/merged_metadata.csv")
metadata<-metadata%>%filter(Run_Tissue_name==slide)
n<-nrow(metadata)
metadata$part<-cut(metadata$fov,
                   #breaks = round(n/30000),#1194A4-PT5,BEME541_Endometriosis,BEME142,BEME543,BEME523-Endometrioma
                   breaks = round(n/10000),
                   labels = F)

# Format and subset
gobj<-createGiottoObject(expression =counts,
                         expression_feat="protein",
                         cell_metadata = metadata)
# 1st network
gobj<-setSpatialLocations(gobj,
                          createSpatLocsObj(metadata%>%select(x_slide_mm,y_slide_mm,cell_ID)),
                          name = 'raw')
gobj = createSpatialGrid(gobject = gobj,
                         sdimx_stepsize = .5,
                         sdimy_stepsize = .5,
                         minimum_padding = 0)

#separate tissue in manageable parts
i<-metadata%>%split(metadata$part)%>%lapply(function(x) x$cell_ID)
gobj<-lapply(i, function(x)  subsetGiotto(gobj,cell_ids = x))
rm(metadata)
gc()

hmrf_folder <- file.path("docs/HMRF_results")
if (!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = TRUE)
for (i in 1:length(gobj)){
  # HMRF requires a fully connected network!
  gobj[[i]] <- createSpatialNetwork(gobj[[i]],
                                    minimum_k = 2)
  #possible spatial genes
  cor_genes <- detectSpatialCorFeats(
    gobj[[i]],expression_values = "raw",
    method = "network",
    spatial_network_name = "Delaunay_network")
  temp<-cor_genes$cor_DT%>%
    select(feat_ID,variable,spat_cor)%>%
    pivot_wider(names_from = variable,values_from = spat_cor)%>%
    column_to_rownames("feat_ID")
  ggsave(filename = paste0(slide,"_",i,"_k.png"),
         plot = fviz_nbclust(temp,FUNcluster = hcut,method = "wss",k.max = 20))
  cor_genes<-cor_genes$cor_DT%>%
    filter(spat_cor>.6&feat_ID!=variable)%>%
    select(feat_ID,variable)%>%unlist()%>%unique
  cor_genes<-cor_genes[cor_genes%in%probes]
  
  HMRF_spatial_genes <- doHMRF(
    gobject = gobj[[i]],
    expression_values = "raw",
    spatial_network_name = "Delaunay_network",
    spatial_genes = cor_genes,
    k =8,#arbitrary high
    tolerance = 1e-8,
    betas = c(13, 5, 1),#start,step,n
    output_folder = file.path(hmrf_folder, paste0(slide,"_",i,"_k8")))
  betas<-list.files(paste0(hmrf_folder,
                           paste0("/",slide,"_",i,"_k8/"),
                           "result.spatial.zscore/k_8/"))
  betas<-betas[str_detect(betas,"hmrf")]
  if(length(betas)==0) stop(paste0("part",i,"failed"))
  gobj[[i]] <- addHMRF(gobject = gobj[[i]],
                       HMRFoutput =HMRF_spatial_genes,
                       k = 8,
                       betas_to_add = 13,
                       hmrf_name = "HMRF")
  
}
metadata<-lapply(gobj,function(x) getCellMetadata(gobject = x,output = "data.table"))

i<-lapply(metadata,function(x) {tab<-x%>%
  select(1,ncol(x)-1,ncol(x))
colnames(tab)[3]<-"domain";
tab})
i<-do.call(rbind,i)
i%>%as.data.frame()%>%data.table::fwrite(paste0("docs/",slide,"_domains.csv"))
gc()

i<-i%>%drop_na()%>%unite("part_domain",part:domain)
i<-i%>%split(i$part_domain)%>%lapply(function(y) y$cell_ID)
domain_signature<-sapply(i,function(x) rowMeans(counts[,x]))
end.time <- Sys.time()
end.time - start.time

write_csv(domain_signature%>%as.data.frame()%>%rownames_to_column("probe"),
          file =paste0("docs/",slide,"_metagenes.csv") )