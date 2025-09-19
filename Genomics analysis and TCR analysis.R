library(Seurat)
dir<-dir()
scRNAlist <- list()
for(i in 1:length(dir)){
counts <- Read10X(data.dir = dir[i])
scRNAlist[[i]] <- CreateSeuratObject(counts, project=dir[i], min.cells=3, min.features = 200)
} 
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
scRNA <- subset(scRNA, subset =nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10) 
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)%>%ScaleData(vars.to.regress=c(cc.genes$s.genes,cc.genes$g2m.genes))%>%RunPCA()
scRNA  <- FindNeighbors(scRNA, dims=1:20) %>% FindClusters(resolution=0.3)
scRNA  <- scRNA %>% RunUMAP(dims=1:20)
##VDJanalysis##
library(scRepertoire)
dir<-dir()
TCR <- list()
for(i in 1:length(dir)){
TCR[[i]] <- read.csv(file=dir[i],stringsAsFactors =F)
}
combined <- combineTCR(TCR, 
                       samples = c("BM-R-25-aCCR5", "BM-R-25-aPD1", "BM-R-25-PD1-aCCR5", "BM-R-25-PDO"))
seurat <- combineExpression(combined, scRNA, 
                            cloneCall="gene", 
                           # group.by = "Sample", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
