#!/path/to/Rscript --vanilla

#cellcols<-c("#BE0032","#66c2a5","cornflowerblue","#f46d43","#8DB600",
#            "#fdae61","#a95f96","hotpink","#0067A5","#604E97","#97604e")
#names(cellcols)<-c("Epithelial","Endothelial cells","Stroma","Myeloid cells","B cells",
#                   "Plasma cells","CD4 T cells","unknown","CD8 T cells","NK cells","Non-immune")


library(tidyverse)
library(Seurat)

matfile<-commandArgs(trailingOnly=TRUE)

counts<-data.table::fread(matfile)
matfile<-str_remove_all(matfile,".+\\/|.csv")
#counts<-data.table::fread("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/PE.csv")
meta<-vroom::vroom(paste0("/media/Lawrenson_Lab_NAS/uthscsa/group_data/CosMx_temp/",matfile,"_meta_bgclus.csv"))#%>%filter(cell_ID%in%counts$cell_ID)
gc()

counts<-counts%>%column_to_rownames("cell_ID")%>%t()
negavg<-colMeans(counts[c("Ms IgG1","Rb IgG"),])
normi<-t(t(counts)-negavg)
normi<-asinh(normi/5)
rm(negavg)
gc()

normi<-sva::ComBat(dat = normi,batch = str_remove(colnames(normi),"c_1.+"),
                   BPPARAM = BiocParallel::MulticoreParam(workers = 2))

normi<-normi%>%t()%>%as.data.frame()%>%rownames_to_column("cell_ID")
normi<-meta%>%select(cell_ID,merged_corrected,cl_type)%>%right_join(normi)%>%drop_na()
ref<-normi%>%group_by(merged_corrected)%>%
  summarise(across(c(CD45,EpCAM,FABP4,CD3,CD38,CD20,CD4,CD8,CD56,CD16,SMA),
                   mean))%>%
  column_to_rownames("merged_corrected")
j<-rownames(counts)
j<-j[str_detect(j,"Chann|IgG",negate = T)]

recluster<-function(index){
  submat<-counts[,index]
  submeta<-meta%>%filter(cell_ID%in%index)%>%
    select(cell_ID,cl_type,merged_corrected)
  submat<-submat[,submeta$cell_ID]
  sc<-CreateSeuratObject(counts = submat,meta.data = submeta)
  sc<-NormalizeData(sc)
  sc<-ScaleData(sc)
  sc<-RunPCA(sc,features = j)
  gc()
  sc<-harmony::RunHarmony(sc,group.by.vars="orig.ident")
  sc<-FindNeighbors(sc,reduction = "harmony")
  resol<-ifelse(ncol(sc)>1e5,5,1)
  sc<-FindClusters(sc,resolution = resol) 
  
  res<-sc@meta.data%>%count(seurat_clusters,merged_corrected)%>%
    group_by(seurat_clusters)%>%slice_max(n)%>%
    mutate(cl_corrected=merged_corrected)%>%select(-merged_corrected)
  res<-sc@meta.data%>%select(cell_ID,seurat_clusters)%>%right_join(res)%>%select(-n)
  return(res)}

#correct epithelia
i<-normi%>%filter(cl_type!="Epithelial"&EpCAM>ref["Epithelial","EpCAM"])%>%
  select(cell_ID)%>%unlist()
reassigned<-recluster(i)
gc()
i<-normi%>%filter(cl_type=="Myeloid cells"&EpCAM>ref["Epithelial","EpCAM"])%>%
  select(cell_ID)%>%unlist()
reassigned<-recluster(i)%>%bind_rows(reassigned)
gc()

#correct stroma
i<-normi%>%filter(cl_type=="Stroma"&
                    CD45>mean(ref[str_detect(rownames(ref),"Epi|Stro",negate = T),"CD45"]))%>%
  select(cell_ID)%>%unlist()
reassigned<-recluster(i)%>%bind_rows(reassigned)
gc()

#correct B cells
i<-normi%>%filter(cl_type=="B cells"&
                    CD3>mean(ref[c("CD8 T cells","CD4 T cells"),"CD3"]))%>%
  select(cell_ID)%>%unlist()
reassigned<-recluster(i)%>%bind_rows(reassigned)
gc()
#correct Plasma
#i<-normi%>%filter(cl_type%in%c("B cells")&
#                    CD38>ref["Plasma cells","CD38"])%>%
#  select(cell_ID)%>%unlist()
#reassigned<-recluster(i)%>%bind_rows(reassigned)
#gc()

#correct CD8 cells
i<-normi%>%filter(cl_type%in%c("CD8 T cells")&
                    CD4>ref["CD4 T cells","CD4"])%>%
  select(cell_ID)%>%unlist()
reassigned<-recluster(i)%>%bind_rows(reassigned)
gc()
#correct CD4 cells
i<-normi%>%filter(cl_type=="CD4 T cells"&
                    (CD8>ref["CD8 T cells","CD8"]|CD16>ref["Myeloid cells","CD16"]))%>%
  select(cell_ID)%>%unlist()
reassigned<-recluster(i)%>%bind_rows(reassigned)
gc()

#correct NK cells
i<-normi%>%filter(cl_type%in%c("Plasma cells")&
                    CD56>ref["NK cells","CD56"])%>%
  select(cell_ID)%>%unlist()
reassigned<-recluster(i)%>%bind_rows(reassigned)
gc()


reassigned<-unique(reassigned)
reassigned<-reassigned%>%add_count(cell_ID)%>%
  mutate(cl_corrected=ifelse(n>1,"unknown",cl_corrected))%>%
  select(-n)%>%distinct()
data.table::fwrite(reassigned,
                   paste0(matfile,"_clcorrected.csv"))