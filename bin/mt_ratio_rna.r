#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

library(dplyr)
library(Seurat)
library(Matrix)

matrix_dir<-as.character(commandArgs(TRUE)[1])
genome<- as.character(commandArgs(TRUE)[2])
summaryFile <- as.character(commandArgs(TRUE)[3])

### Ratio GeneVSumi & %MT for mqc plots
#--------------------------------------
sparseMtx <- Read10X(data.dir = matrix_dir)
resume<-read.table(summaryFile, sep=',', row.names=TRUE, col.names=FALSE)

# Ratio
Ratio<-cbind(rownames(resume), resume)
colnames(Ratio)<-c("Samples", "Number of genes", "Number of UMIs")
write.table(Ratio, "RatioPerCell.csv", sep=',', row.names=FALSE, col.names=FALSE)

# MT
umiMatrix <- CreateSeuratObject(counts = sparseMtx, min.features = 0)
if(genome=="hg38" || genome=="hg19"){
    umiMatrix[["percent.mt"]] <- PercentageFeatureSet(umiMatrix, pattern = "^MT-")
}
# if mouse genome
if(genome=="mm10" || genome=="mm9" || genome=="dmelr6.28"){
    umiMatrix[["percent.mt"]] <- PercentageFeatureSet(umiMatrix, pattern = "^mt-")
}

MT<-cbind(rownames(umiMatrix[["percent.mt"]]), Ratio$`Number of genes`, umiMatrix[["percent.mt"]])
colnames(MT)<-c("Samples", "Number of genes", "% Mitochondrial genes")
write.table(MT, "MtGenePerCell.csv", sep=',', row.names=FALSE, col.names=FALSE)


