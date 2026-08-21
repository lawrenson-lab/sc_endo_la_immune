#!/path/to/Rscript --vanilla

library(tidyverse)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
#reticulate::use_condaenv("/home/msoledad/.virtualenvs/renv/bin/python3")

matfile<-commandArgs(trailingOnly=TRUE)

counts<-data.table::fread(matfile)
matfile<-str_remove_all(matfile,".+\\/|.csv")
counts<-counts%>%column_to_rownames("cell_ID")%>%t()

negavg<-colMeans(counts[c("Ms IgG1","Rb IgG"),])
normi<-t(t(counts)-negavg)
normi<-asinh(normi/5)
nobatch<-sva::ComBat(dat = normi,batch = str_remove(colnames(normi),"c_1.+"),
                     BPPARAM = BiocParallel::MulticoreParam(workers = 2))
#normi<-nobatch/apply(nobatch,1,max)
#
sctype_annotate<-function(mrkrs,mtrx){
  gs_positive<-mrkrs%>%
    split(mrkrs$Annotation)%>%
    lapply(function(x) x$Positive)
  gs_negative<-lapply(gs_positive,function(x) 
    unique(unlist(gs_positive))[!unique(unlist(gs_positive))%in%x])
  cell_score <- sctype_score(mtrx, 
                             scaled = F, 
                             gs = gs_positive, 
                             gs2 = gs_negative,
                             gene_names_to_uppercase = F)
  res<-as.data.frame(cbind(cell_ID=colnames(cell_score),
                           cell_type=rownames(cell_score)[apply(cell_score,2,which.max)],
                           score=apply(cell_score,2,max),
                           t(cell_score)))
  return(res)}
markers<-read_tsv("cosmx_markers.tsv")
markers<-markers%>%select(-Negative)%>%
  separate_rows(Positive,sep=',')%>%
  mutate(Positive=str_remove(Positive,' '))
main<-markers%>%filter(Annotation%in%c("Stroma","Immune","Epithelial","Adipose"))
markers<-markers%>%filter(Annotation!="Immune")

assignment<-sctype_annotate(main,nobatch)
assignment<-assignment%>%mutate(across(score:Stroma,as.numeric))
assignment<-assignment%>%mutate(cell_type=ifelse(score<0,"unknown",cell_type))
colnames(assignment)[-1]<-paste0("main_",colnames(assignment)[-1])
res<-sctype_annotate(markers,nobatch)
res<-res%>%mutate(across(score:Stroma,as.numeric))
assignment<-assignment%>%left_join(res)

assignment<-assignment%>%
  mutate(merged_type=ifelse(main_cell_type%in%c("Immune","unknown"),cell_type,
                            main_cell_type))#some Immune

assignment<-assignment%>%
  mutate(merged_corrected=ifelse(str_detect(merged_type,"Adip|Stro|Endot"),"Stroma",
                                 merged_type))

assignment%>%data.table::fwrite(paste0(matfile,"_sctype.csv"))
