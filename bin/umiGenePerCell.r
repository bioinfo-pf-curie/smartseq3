#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

library(reshape2)

#dir_matrices<-as.character(commandArgs(TRUE)[1])

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

if(length(listFile)==1){
    #If only one cell
    nbUMIs<-sum(matrixFinal[,2])
    nbGenes<-nrow(matrixFinal)
}else{
    nbUMIs<-colSums(matrixFinal[,-1])
    nbGenes<-apply(matrixFinal[,-1], 2, function(x) length(which(x>0)))
}

lg_nbUMIs<-melt(as.matrix(nbUMIs))
lg_nbUMIs<-lg_nbUMIs[,-2]
colnames(lg_nbUMIs)<-c("Cells", "Number of umis")
write.table(lg_nbUMIs, "nbUMIPerCell.csv",
            sep=',', row.names=FALSE, col.names=TRUE)

lg_nbGenes<-melt(as.matrix(nbGenes))
lg_nbGenes<-lg_nbGenes[,-2]
colnames(lg_nbGenes)<-c("Cells", "Number of genes")
write.table(lg_nbGenes, "nbGenePerCell.csv",
            sep=',', row.names=FALSE, col.names=TRUE)
