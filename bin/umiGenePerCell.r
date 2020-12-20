#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))


umiMatrix<- as.character(commandArgs(TRUE)[1])
prefix = as.character(commandArgs(TRUE)[2])

matrix<-read.table(umiMatrix, header=TRUE)

nbGenes<-nrow(matrix)
write.table(nb_Genes, paste0(as.character(prefix), "_nbGenePerCell_mqc.csv"),
       sep=',', row.names=FALSE, col.names=FALSE)

nbUmis<-sum(matrix[2])
write.table(nbUmis, paste0(as.character(prefix), "_nbUMIPerCell_mqc.csv"),
            sep=',', row.names=FALSE, col.names=FALSE)
