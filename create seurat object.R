library(Seurat)
dir<-dir()
scRNAlist <- list()
for(i in 1:length(dir)){
counts <- Read10X(data.dir = dir[i],gene.column = 1)
scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample[i], min.cells=3, min.features = 200)
} 
##qdoublet filter##
library(DoubletFinder)
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:50, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(0.05*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:50, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  return(data)
}
for(i in 1:length(prostate)){
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset =nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20) 
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]])%>%ScaleData()%>%RunPCA()
  scRNAlist[[i]] <- FindNeighbors(scRNAlist[[i]], dims=1:30) %>% FindClusters(resolution=0.2)
  scRNAlist[[i]] <- scRNAlist[[i]] %>% RunTSNE(dims=1:30,check_duplicates=F) %>% RunUMAP(dims=1:30)
  scRNAlist[[i]] <-Find_doublet(scRNAlist[[i]])   
  c <- grep("pANN_",colnames(scRNAlist[[i]]@meta.data))
  scRNAlist[[i]]@meta.data <- scRNAlist[[i]]@meta.data[,-c]
  scRNAlist[[i]]<-subset(scRNAlist[[i]],subset=doublet_info=="Singlet")
}
##analysis as a whole##
scRNA<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
HSP_genes <- rownames(scRNA)[grep("^HSP",rownames(scRNA))]
RLS_genes <- rownames(scRNA)[grep("^RP[SL]",rownames(scRNA))]
mt.genes_1 <- rownames(scRNA)[grep("^MT-",rownames(scRNA))]
mt.genes <- rownames(scRNA)[grep("^MT1",rownames(scRNA))]
filter.gene<-unique(c(HSP_genes,mt.genes,cc.genes$g2m.genes,cc.genes$s.genes))
scRNA <-scRNA%>% ScaleData(vars.to.regress=filter.gene)   
features<-VariableFeatures(scRNA)
features<-setdiff(features,filter.gene)
scRNA <-scRNA%>%RunPCA(features=features)
p=ElbowPlot(scRNA,ndims=50)
n=ndims
scRNA <- FindNeighbors(scRNA, dims=1:n) %>% FindClusters(resolution=0.1)
scRNA <-RunUMAP(scRNA,dims=1:n)
##analysis by cancer type##
scRNAlist<-lapply(unique(scRNA$cancertype),function(x) subset(scRNA,subset=cancertype==x))
scRNAlist[[i]]<-{
HSP_genes <- rownames(scRNAlist[[i]])[grep("^HSP",rownames(scRNAlist[[i]]))]
RLS_genes <- rownames(scRNAlist[[i]])[grep("^RP[SL]",rownames(scRNAlist[[i]]))]
mt.genes_1 <- rownames(scRNAlist[[i]])[grep("^MT-",rownames(scRNAlist[[i]]))]
mt.genes <- rownames(scRNAlist[[i]])[grep("^MT1",rownames(scRNAlist[[i]]))]
filter.gene<-unique(c(HSP_genes,mt.genes,cc.genes$g2m.genes,cc.genes$s.genes))
scRNAlist[[i]] <-scRNAlist[[i]]%>% ScaleData(vars.to.regress=filter.gene)   
features<-VariableFeatures(scRNAlist[[i]])
features<-setdiff(features,filter.gene)
scRNAlist[[i]] <-scRNAlist[[i]]%>%RunPCA(features=features)
}
## determine the number of dims##
plot<-list()
for (i in 1:length(unqiue(scRNA$cancertype)))
{
  plot[[i]]=ElbowPlot(scRNAlist[[i]],ndims=50)
  temp<-plot[[i]]
ggsave(paste0("BM_MGI_elbowplot_",unqiue(scRNA$cancertype)[[i]],".png"),temp,width=5,height=5)
}
n=dims
scRNAlist[[i]]<-{
scRNAlist <- FindNeighbors(scRNAlist, dims=1:n) %>% FindClusters(resolution=0.1)
scRNAlist <-RunUMAP(scRNAlist,dims=1:n)}