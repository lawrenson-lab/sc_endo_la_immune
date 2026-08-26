#!/usr/bin/env Rscript

######################### load stuff###############################
library(tiledb,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/")
library(tiledbsc,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/")
library(tidyverse,lib.loc="/home/msoledad/cosmx/renv/library/R-4.3/x86_64-pc-linux-gnu/")

file<-commandArgs(trailingOnly=TRUE)
tiledb_scdataset <- SOMACollection$new(uri = file,
                                       verbose = F)#F or it'll print a ton of filenames
metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()

slide<-metadata$Run_Tissue_name[1]

counts <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix()
counts<-counts[!str_detect(rownames(counts),"-DNA|-G|-Membrane"),]

normi<-as.matrix(asinh(counts/5))

normi%>%as.data.frame%>%
  rownames_to_column("probe")%>%
  data.table::fwrite(paste0("docs/",slide,"asinh.csv"))