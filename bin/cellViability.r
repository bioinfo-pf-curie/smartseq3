#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

library(plotrix)
library(plyr)
library(dplyr)
library(Seurat)
library(tibble)
library(reshape2)
library(ggplot2)

dir_matrices<-as.character(commandArgs(TRUE)[1])
dir_res_10X<-as.character(commandArgs(TRUE)[2])

##### Merge Matrices
####----------------------------------------

listFile<-list.files(path = dir_matrices)

i=1
for (file in listFile ){
    matrix<-read.table(file = paste(dir_matrices, file, sep = ""), header=TRUE)
    #matrix<-read.table(file = paste("matrices/mouse/", file, sep = ""), header=TRUE)
    
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
# Transform into long format
longMatx<-melt(matrixFinal)
longMatx<-longMatx[which(longMatx$value>0),]

matrixFinal<-column_to_rownames(matrixFinal, var = "gene")

# Normalize each count by the total number of umis per cell *10^5
#normData<-data.frame(NormalizeData(matrixFinal, normalization.method="RC" , scale.factor = 100000))
# log transformation
#normLogData<-log10(normData + 1)
# transform into long format
#normLogDataLong<-melt(as.matrix(normLogData))
#normLogDataLong<-normLogDataLong[which(normLogDataLong$value>0),]

# Make resume data 
resume<-data.frame(nb_UMIs = colSums(matrixFinal))
resume$nb_Genes<-apply(matrixFinal, 2,  function(x) length(which(x>0)))

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


#--------------------------
# Nb UMI per Gene
#=> par cellules (cumulative dist) mais pas de zéro (col1=sample; col2 = counts; col3=gene name )

#nbUMIs_perGene<-dcast(longMatx, formula = longMatx2$value~longMatx2$variable, value.var = "value", fun.aggregate = length)

# sans log10
#longMatx2<-select(longMatx, -gene)
#nbUMIs_perGene<-dcast(longMatx2, formula = longMatx2$value~longMatx2$variable, value.var = "value", fun.aggregate = length)
#names(nbUMIs_perGene)[[1]]<-"x axis"
#nbUMIs_perGene<-nbUMIs_perGene[,-1]
# si test avec log10
#nbUMIs_perGene$log10<-log10(nbUMIs_perGene$`longMatx2$value`)

# si test avec 
#nbUMIs_perGene_max50<-nbUMIs_perGene[c(1:100),]

#write.csv(nbUMIs_perGene, "HistUMIperGene.mqc", row.names = FALSE)

#write.table(t_nbUMIs_perGene, "HistUMIperGene.mqc", sep=',', row.names=TRUE, col.names=FALSE)

#--------------------------
# Nb UMI & Gene per cell
# Jitter sur le même graphe avec médiane (image fixe)

# Les 2 sur le même mais en log10
long_resume<-melt(as.matrix(resume))
colnames(long_resume)<-c("Samples", "Var2", "value")
jpeg(file="jitter_nbUMI_nbGenes_mqc.jpeg")
ggplot(data = long_resume, aes(x = Var2, y = value)) + theme_bw() +
        geom_jitter(mapping = aes(colour = Samples), width = .1)  +
        stat_summary(fun=median, geom="point", shape=18,  size=3) +
        xlab("") + ylab("Counts (log10)") + 
        scale_y_log10()
dev.off()
#ggsave("jitter_nbUMI_nbGenes.tiff", units="in", width=5, height=4, dpi=300)

#---------------

### ratio GeneVSumi & %MT
#------------
#=> scatter en renommant les axes:

Ratio<-cbind(rownames(resume), resume)
colnames(Ratio)<-c("Samples", "Number of genes", "Number of UMIs")
#write.csv(Ratio, "RatioPerCell_mqc.csv", row.names = FALSE )
# => donne un graphe still heatmap 

write.table(Ratio, "RatioPerCell_mqc.csv", sep=',', row.names=FALSE, col.names=FALSE)


umiMatrix <- CreateSeuratObject(counts = sparseMtx, min.features = 0)
umiMatrix[["percent.mt"]] <- PercentageFeatureSet(umiMatrix, pattern = "^MT-")
MT<-cbind(rownames(umiMatrix[["percent.mt"]]), Ratio$`Number of genes`, umiMatrix[["percent.mt"]])
colnames(MT)<-c("Samples", "Number of genes", "% Mitochondrial genes")
#write.csv(MT, "MtGenePerCell_mqc.csv", row.names = FALSE )

write.table(MT, "MtGenePerCell_mqc.csv", sep=',', row.names=FALSE, col.names=FALSE)


