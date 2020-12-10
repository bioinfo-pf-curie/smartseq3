#!/usr/bin/env Rscript

library(plotrix)
library(plyr)
library(dplyr)
library(Seurat)
library(tibble)
library(reshape2)
library(ggplot2)
.libPaths()
.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))
.libPaths()
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
        
        sample<-strsplit(x = file ,split = "ID") [[1]][1]
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
normData<-data.frame(NormalizeData(matrixFinal, normalization.method="RC" , scale.factor = 100000))
# log transformation
normLogData<-log10(normData + 1)
# transform into long format
normLogDataLong<-melt(as.matrix(normLogData))
normLogDataLong<-normLogDataLong[which(normLogDataLong$value>0),]

# Make resume data 
samples<-colnames(matrixFinal)
genes<-vector()
for(samp in samples){
    genes<-append(genes, nrow(longMatx[which(longMatx$variable==samp),]) )
}

resume<-data.frame(nb_UMIs = colSums(matrixFinal), NormLog_nb_UMIs = round(colSums(normLogData)) , nb_Genes=genes)
write.csv(resume, "resume.csv", row.names = T)


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

### Counts
#------------
hist_nbUMIperGene<-hist(normLogDataLong$value, xlab = "# UMIs (Log10Normalized)", ylab = "# Genes", main = "Number of UMIs per genes")
hist_nbUMIperCell<-hist(resume$NormLog_nb_UMIs, xlab = "# UMIs (Log10Normalized)", ylab = "# Cell")
hist_nbGenesPerCell<-hist(resume$nb_Genes, xlab = "# Genes", ylab = "# Cell")
wh_UMI<-weighted.hist(resume$NormLog_nb_UMIs,
              w=resume$NormLog_nb_UMIs, 
              main="Weighted distribution of umis per barcodes",
              xlab="#UMIs per cell (log10)",
              ylab="#UMIs", breaks = 5)

create_df<-function(list){
    df<-data.frame(list$breaks)
    df<-df[-1,]
    df<-cbind(df, list$counts)
    colnames(df)<-c("breaks.x", "counts.y")
    return(df)
}

nbUMIperGene<-create_df(hist_nbUMIperGene)
nbUMIperCell<-create_df(hist_nbUMIperCell)
nbGenesPerCell<-create_df(hist_nbGenesPerCell)
whUMI<-create_df(wh_UMI)

write.csv(nbUMIperGene, "HistUMIperGene.csv")
write.csv(nbUMIperCell, "HistUMIperCell.csv")
write.csv(nbGenesPerCell, "HistGenesPerCell.csv")
write.csv(wh_UMI, "weightedHistUMI.csv")

### ratio GeneVSumi & %MT
#------------

umiMatrix <- CreateSeuratObject(counts = sparseMtx, project = "smartSeq3", min.features = 0)
umiMatrix[["percent.mt"]] <- PercentageFeatureSet(umiMatrix, pattern = "^MT-")

plotRatio<-FeatureScatter(umiMatrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Ratio<-data.frame(plotRatio[["data"]][["nFeature_RNA"]])
Ratio$umi<-plotRatio[["data"]][["nCount_RNA"]]
colnames(Ratio)<-c("#Genes", "#UMIs")
write.csv(Ratio, "UmiGenePerCell.csv")

plotMT<-FeatureScatter(umiMatrix, feature1 = "nFeature_RNA", feature2 = "percent.mt")
MT<-data.frame(plotMT[["data"]][["nFeature_RNA"]])
MT$mt<-plotMT[["data"]][["percent.mt"]]
colnames(MT)<-c("#Genes", "%MT")
write.csv(MT, "MtGenePerCell.csv")


