#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

library(reshape2)

dir_matrices<-as.character(commandArgs(TRUE)[1])

listFile<-list.files(path = dir_matrices)

#dir_matrices<-"matrices/"
i=1
for (file in listFile ){
    matrix<-read.table(file = paste(dir_matrices, file, sep = ""), header=TRUE)
    if( nrow(matrix)!=0 ){
        # If repeated gene names (different chr (often Y_RNA))
        matrix<-aggregate(matrix$count ~ matrix$gene, FUN=sum)
        
        sample<-strsplit(x = file ,split = "_") [[1]][1]
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

nbUMIs<-colSums(matrixFinal[,-1])
lg_nbUMIs<-melt(as.matrix(nbUMIs))
lg_nbUMIs<-lg_nbUMIs[,-2]
colnames(lg_nbUMIs)<-c("Cells", "Number of umis")
write.table(lg_nbUMIs, "nbUMIPerCell.csv",
           sep=',', row.names=FALSE, col.names=TRUE)


nbGenes<-apply(matrixFinal[,-1], 2, function(x) length(which(x>0)))
lg_nbGenes<-melt(as.matrix(nbGenes))
lg_nbGenes<-lg_nbGenes[,-2]
colnames(lg_nbGenes)<-c("Cells", "Number of genes")
write.table(lg_nbGenes, "nbGenePerCell.csv",
            sep=',', row.names=FALSE, col.names=TRUE)


#################
# si matrice une à une:

# umiMatrix<- as.character(commandArgs(TRUE)[1])
# prefix = as.character(commandArgs(TRUE)[2])
#matrix<-read.table(umiMatrix, header=TRUE)

#nbGenes<-nrow(matrix)
#colnames(nbGenes)<-c("Cell", "Number of genes")
#write.table(nbGenes, paste0(as.character(prefix), "nbGenePerCell_mqc.csv", collapse = "_"),
#            sep=',', row.names=FALSE, col.names=TRUE)


# nbUmis<-sum(matrix[2])
# colnames(nbUmis)<-c("Cell", "Number of UMIs")
# write.table(nbUmis, paste0(as.character(prefix), "nbUMIPerCell_mqc.csv", collapse = "_"),
#              sep=',', row.names=FALSE, col.names=TRUE)

#countUMIGene<-cbind(prefix, nbGenes, nbUmis)
#colnames(countUMIGene)<-c("Cell", "Number of genes", "Number of UMIs")
#
# write.table(countUMIGene, paste0(as.character(prefix), "_countPerCell_mqc.csv"),
#             sep=',', row.names=FALSE, col.names=FALSE)
# 
