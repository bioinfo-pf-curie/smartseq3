#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

library(reshape2)

umiMatrix<- as.character(commandArgs(TRUE)[1])
prefix = as.character(commandArgs(TRUE)[2])

matrix<-read.table(umiMatrix, header=TRUE)

# TEST:
# setwd("~/Documents/Curie/SC_platform/Analysis/scRNAseq/1.Demultiplexing/smartSeq3/tests/create10X/")
# matx<-read.table(file =  "matrices/toy_L384T27_umi_Counts.tsv.gz", header = T)
# matrix<-matx
# prefix<-"toy_L384T27_umi"

nbUMIs_perGene<-dcast(matrix, formula = matrix$count ~ prefix, value.var = "count", fun.aggregate = length)
colnames(nbUMIs_perGene)<-c("x","y")

write.table(nbUMIs_perGene, paste0(as.character(prefix), "_HistUMIperGene_mqc.csv"),
            sep=',', row.names=FALSE, col.names=FALSE)