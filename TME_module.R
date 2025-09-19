library(Seurat)
library(ggplot2)
library(ggpubr)
##TME module analysis##
BM<-readRDS("BM_MGI_RNA.rds")
BM_1<-subset(BM,subset=celltype!="Tumor_cell")
cellratio<-prop.table(table(BM_1$cellcluster,BM_1$orig.ident),margin=2)
cor<-cor(cellratio,method="pearson")
color = colorRampPalette(colors = c("Navy","white","#8F2626"))(100)
pdf("BM_BMt_cancertype_cellcluster_cor_integrated_total.pdf",width=12,height=12)
corrplot(cor,order='hclust',method='color',hclust.method="average",col=color,tl.cex=0.7)
##distribution of TME module##
cellratio<-prop.table(table(BM_1$module,BM_1$orig.ident),margin=2)
write.csv(cellratio,"BM_BMT_module_cellratio_orig.ident.csv")
cellratio<-read.csv("BM_BMT_module_cellratio_orig.ident.csv",header=T,row.names=1)
cellratio<-as.data.frame(t(cellratio))
meta<-BM_1@meta.data[c("orig.ident","cancertype")]
meta<-meta[!duplicated(meta$orig.ident),]
rownames(meta)<-meta$orig.ident
meta<-meta[rownames(cellratio),]
cellratio$cancertype<-meta$cancertype
ratio<-lapply(unique(cellratio$cancertype),function(x) subset(cellratio,subset=cancertype==x))
mean_ratio<-lapply(ratio,function(x) apply(x[,-ncol(x)],2,mean))
mean_ratio<-do.call("cbind",mean_ratio)
colnames(mean_ratio)<-unique(cellratio$cancertype)