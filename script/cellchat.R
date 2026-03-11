library(CellChat)
library(tidyverse)
library(Seurat)

endo_cluster <- readRDS("data/fixed_annotated_aux.seurat.rds")
endo_cluster<-subset(endo_cluster,cells=which(endo_cluster$Class=="Peritoneal endometriosis"))
#Making group_id(s) for cell chat comparison 


run_subset_cell_chat <- function(cell_chat_data,group_id){
  meta_table <- cell_chat_data@meta.data
  cell_chat_data <- GetAssayData(cell_chat_data, assay = "RNA", slot = "data")
  cell_chat_data <- createCellChat(object = cell_chat_data, meta = meta_table, group.by = group_id)
  
  print("Cells start Chatting !")
  
  cell_chat_data@DB <- CellChatDB.human
  
  #pre-processing
  print("pre-processing is now commenced")
  cell_chat_data <- cell_chat_data %>%
    subsetData %>%
    identifyOverExpressedGenes %>%
    identifyOverExpressedInteractions #%>% 
  #projectData(PPI.human) 
  
  print("pre-processsing done -let's start processing!")
  #post-processing
  
  cell_chat_data <- cell_chat_data %>%  
    computeCommunProb(.,raw.use = TRUE,population.size = TRUE) %>%
    filterCommunication(.,min.cells = 10) %>%
    computeCommunProbPathway %>%
    aggregateNet
  
  
}

make_n_merge_cell_chat <- function(seurat_obj,subset_var,group_id) {
  seurat_obj <- SplitObject(seurat_obj,split.by =subset_var)
  object.list <- lapply(seurat_obj,run_subset_cell_chat,group_id)
  cellchat <- mergeCellChat(object.list, add.names = names(seurat_obj))
  cellchat <- c(cellchat, object.list)                         
}

endo_cluster$has_LA<-endo_cluster$LA_per_Sample>0
endo_cluster<-subset(endo_cluster,
                     cells=which(endo_cluster$harmonized_major_plus_immune_cell_type != "Exclude"))
endo_cluster@meta.data<-endo_cluster@meta.data%>%
  mutate(nosub=str_remove(harmonized_major_plus_immune_cell_type," \\(.+"))%>%
  mutate(nosub=str_replace(nosub,"Memmo","Memo")%>%str_remove("Early "))

cell_chat_data <- make_n_merge_cell_chat(endo_cluster,"has_LA", "nosub")
saveRDS(cell_chat_data,"data/CellChat_PE_hasLA_exclude.rds")

