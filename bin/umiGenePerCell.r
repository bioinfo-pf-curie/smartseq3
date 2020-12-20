#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))


umiMatrix<- as.character(commandArgs(TRUE)[1])
prefix = as.character(commandArgs(TRUE)[2])

matrix<-read.table(umiMatrix, header=TRUE)

nbGenes<-nrow(matrix)
# write.table(nbGenes, paste0(as.character(prefix), "_nbGenePerCell_mqc.csv"),
#        sep=',', row.names=FALSE, col.names=FALSE)

nbUmis<-sum(matrix[2])

countUMIGene<-cbind(prefix, nbGenes, nbUmis)

colnames(countUMIGene)<-c("Cell", "Number of genes", "Number of UMIs")

write.table(countUMIGene, paste0(as.character(prefix), "_countPerCell_mqc.csv"),
            sep=',', row.names=FALSE, col.names=TRUE)

