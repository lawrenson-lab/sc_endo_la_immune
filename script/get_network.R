#!/usr/bin/env Rscript
library(tidyverse)
library(geometry)

sampl<-commandArgs(trailingOnly=TRUE)

metadata<-data.table::fread("docs/endometrium_meta_bgclus.csv")#PE_meta_bgclus.csv
metadata<-metadata%>%
  filter(Run_Tissue_name==sampl)

spatial_locations<-metadata%>%distinct(cell_id,x_slide_mm,y_slide_mm)
cell_ID<-spatial_locations%>%select(cell_id)%>%
  mutate(index=1:nrow(spatial_locations))%>%deframe()
spatial_locations<-spatial_locations%>%column_to_rownames("cell_id")
delaunay_triangles  <- geometry::delaunayn(p = spatial_locations,options = "Pp")
delaunay_edges <- rbind(delaunay_triangles[ ,c(1,2)],
                        delaunay_triangles[ ,c(1,3)],
                        delaunay_triangles[ ,c(2,3)])%>%
  unique()%>%
  as.data.frame()%>%
  setNames(c("from","to"))

delaunay_DT = data_frame(from = names(cell_ID[delaunay_edges$from]),
                         to = names(cell_ID[delaunay_edges$to]),
                         xstart = spatial_locations[delaunay_edges$from, "x_slide_mm"],
                         ystart = spatial_locations[delaunay_edges$from, "y_slide_mm"],
                         xend = spatial_locations[delaunay_edges$to, "x_slide_mm"],
                         yend = spatial_locations[delaunay_edges$to, "y_slide_mm"])
delaunay_DT<-delaunay_DT%>%mutate(xlen=xstart-xend,
                                  ylen=yend-ystart,
                                  distance=sqrt(xlen**2+ylen**2))
delaunay_DT%>%select(from,to,distance)%>%write_csv(file = paste0("docs/",sampl,"delaunay_0.1.csv"))
