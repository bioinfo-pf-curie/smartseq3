#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

library(plotrix)
library(plyr)
library(dplyr)
library(Seurat)
library(tibble)
library(reshape2)
library(ggplot2)

dir_res_10X<-as.character(commandArgs(TRUE)[1])

##### Merge Matrices
####----------------------------------------

listFile<-list.files(path = ".", pattern= "*.tsv.gz")

i=1
for (file in listFile ){
    matrix<-read.table(file, header=TRUE)
    if( nrow(matrix)!=0 ){
        # If repeated gene names (different chr (often Y_RNA))
        matrix<-aggregate(matrix$count ~ matrix$gene, FUN=sum)
        
        sample<-strsplit(x = file ,split = "_umi") [[1]][1]
        colnames(matrix)<-c("gene", sample)
        if(i==1){
            matrixFinal<-matrix
        }else{
            matrixFinal<-merge(matrixFinal, matrix, by = "gene", all= TRUE)
        }
        i=2
    }
}

# Replace NA by 0 (row= genes/columns=sample)
matrixFinal[is.na(matrixFinal)]<-0
# Transform into long format
longMatx<-melt(matrixFinal)
longMatx<-longMatx[which(longMatx$value>0),]

matrixFinal<-column_to_rownames(matrixFinal, var = "gene")

# Make resume data 
resume<-data.frame(nb_UMIs = colSums(matrixFinal))
resume$nb_Genes<-apply(matrixFinal, 2,  function(x) length(which(x>0)))
colnames(resume)<-c("UMI_counts", "Gene_counts")

write.table(resume, "resume.txt", sep=',', row.names=TRUE, col.names=FALSE)

#####  Long to sparse matrix:
####-----------------------------------------
library(Matrix)
longMatx$gene=as.factor(longMatx$gene)
longMatx$variable=as.factor(longMatx$variable)
longMatx$value=as.factor(longMatx$value)
sparseMtx <- with(longMatx, sparseMatrix(i=as.numeric(gene),
                                 j=as.numeric(variable),
                                 x=as.numeric(value),
                                 dimnames=list(levels(gene), levels(variable))
)
)

#### Save results like 10X
####-----------------------------------------
library("DropletUtils")
cell.ids <- levels(longMatx[,2])
ngenes <- nrow(sparseMtx)
gene.ids <- paste0("ID", seq_len(ngenes))
gene.symb <-levels(longMatx[,1])

write10xCounts(path = dir_res_10X, sparseMtx, gene.id=gene.ids,
               gene.symbol=gene.symb, barcodes=cell.ids)


#### Create Plots
####-----------------------------------------

### ratio GeneVSumi & %MT
#------------
#=> scatters

Ratio<-cbind(rownames(resume), resume)
colnames(Ratio)<-c("Samples", "Number of genes", "Number of UMIs")
write.table(Ratio, "RatioPerCell.csv", sep=',', row.names=FALSE, col.names=FALSE)


umiMatrix <- CreateSeuratObject(counts = sparseMtx, min.features = 0)
umiMatrix[["percent.mt"]] <- PercentageFeatureSet(umiMatrix, pattern = "^MT-")
MT<-cbind(rownames(umiMatrix[["percent.mt"]]), Ratio$`Number of genes`, umiMatrix[["percent.mt"]])
colnames(MT)<-c("Samples", "Number of genes", "% Mitochondrial genes")
write.table(MT, "MtGenePerCell.csv", sep=',', row.names=FALSE, col.names=FALSE)


