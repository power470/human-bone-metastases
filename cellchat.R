library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(CellChat)
library(NMF)
library(ggalluvial)
scRNA<-readRDS("BM_MGI_RNA.rds")
scRNAlist<-SplitObject(scRNA,split.by="cancertype")
cellchat <- list()
for (i in (1:length(scRNAlist)))
{Idents(scRNAlist[[i]])<-scRNAlist[[i]]$celltype
  DefaultAssay(scRNAlist[[i]])<-"RNA"
cellchat[[i]]<-createCellChat(object=scRNAlist[[i]])
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat[[i]]@DB <- CellChatDB.use
cellchat[[i]] <- subsetData(cellchat[[i]])
cellchat[[i]] <- identifyOverExpressedGenes(cellchat[[i]])
cellchat[[i]] <- identifyOverExpressedInteractions(cellchat[[i]])
cellchat[[i]] <- projectData(cellchat[[i]], PPI.human)
cellchat[[i]] <- computeCommunProb(cellchat[[i]])
cellchat[[i]] <- filterCommunication(cellchat[[i]], min.cells = 20)
cellchat[[i]] <- computeCommunProbPathway(cellchat[[i]])
cellchat[[i]] <- aggregateNet(cellchat[[i]])
groupSize<-as.numeric(table(cellchat[[i]]@idents))}
##LR analysis##
df<-lapply(cellchat,function(x) subsetCommunication(x))
for(i in 1:length(df)){
df[[i]]$cancertype<-unique(scRNA$cancertype)[i]
}
##receptor##
cancertype<-lapply(df,function(x) subset(x,subset=target=="Tumor_cell"))
receptor<-lapply(cancertype,function(x) unique(x[,"receptor"]))
receptor_common<-Reduce(intersect,receptor)
cancertype_common_receptor<-lapply(cancertype,function(x) subset(x,subset=receptor%in%receptor_common))
cancertype_common_receptor<-lapply(receptor_common,function(x) lapply(cancertype_common_receptor,function(y) subset(y,subset=receptor==x)))
interact_mean<-lapply(cancertype_common_receptor,function(x) lapply(x,function(y) mean(y[,"prob"])))
receptor<-interact_mean
receptor<-lapply(receptor, function(x) do.call("rbind",x))
##different receptor##
receptor_all<-lapply(cancertype,function(x) unique(x[,"receptor"]))
receptor_diff<-lapply(receptor_all,function(x) setdiff(x,receptor_common))
receptor_diff_cancertype<-list()
for(i in 1:length(cancertype)){
receptor_diff_cancertype[[i]]<-subset(cancertype[[i]],subset=receptor%in%receptor_diff[[i]])}
receptor_diff_receptor<-lapply(receptor_diff_cancertype,function(x) lapply(unique(x[,"receptor"]),function(y) subset(x,subset=receptor==y)))
interact_mean<-lapply(receptor_diff_receptor,function(x) lapply(x,function(y) mean(y[,"prob"])))
receptor<-interact_mean
receptor<-lapply(receptor, function(x) do.call("rbind",x))
for(i in 1:length(receptor)){
receptor[[i]]<-as.data.frame(receptor[[i]])
receptor[[i]]$cancertype<-unique(cellchat$cancertype)[i]
rownames(receptor[[i]])<-receptor_diff[[i]]
receptor[[i]]$receptor<-rownames(receptor[[i]])
}
receptor<-do.call("rbind",receptor)