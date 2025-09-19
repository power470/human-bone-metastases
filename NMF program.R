library(reshape2)
library(NMF)
library(ggplot2)
library(scales)
##filter gene##
mean_gene<-function(x)
{x<-x@assays$RNA@counts
x<-as.data.frame(x)
x$mean<-apply(x,1,mean)
x<-x[order(x$mean,decreasing=T),]
x<-x[1:7000,-ncol(x)]
x<-scale(x)
x<-data.frame(x)
}
PT2<-mean_gene(BM)
##NMF analysis##
NMF<-list()
ranknumber<-c(4:9)
for (i in 1:length(ranknumber))
{NMF[[i]]<-nmf(PT2,rank=ranknumber[[i]], method="snmf/r", nrun = 10)
}
##nmf programs for each cancer type##
NMF<-list()
for(i in 1:length(COAD)){
NMF[[i]]<- readRDS(file = COAD[i])
} 
nmf_programs<-list()
for (i in 1:length(NMF)){
nmf_programs[[i]]<-lapply(NMF[[i]],function(x){
f <- extractFeatures(x, 50L)
f <- lapply(f, function(y) rownames(x)[y])
do.call("cbind", f)})
}
for (i in 1:length(nmf_programs))
{nmf_programs[[i]]<-do.call("cbind",nmf_programs[[i]])}
##robust NMF programs for each cancer type COAD as a example##
intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) {nmf_programs[[x]][,intra_intersect_max[[x]]>=35]}) 
  names(nmf_sel) <- names(nmf_programs)
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y))))
  inter_filter<-T  
final_filter <- NULL 
sample<-names(nmf_sel)
filter<-function(sample){
for(i in sample) {
    a <- inter_intersect[grep(sample, colnames(inter_intersect), invert = T),grep(sample, colnames(inter_intersect))]
    b <- sort(apply(a, 2, mean), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
    if(inter_filter==T) b <- b[b>=10] # selects programs with a maximum intersection of at least 10
    if(length(b) > 1) {
      c <- names(b[1]) 
      for(y in 2:length(b)) {
        if(max(inter_intersect[c,names(b[y])]) <=25) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
 return(final_filter)    
  }}
final<-lapply(sample,function(x) filter(x))
final<-unlist(final)
nmf_sel_unlist <- nmf_sel_unlist[,final]
write.csv(nmf_sel_unlist,"BM_BMt_COAD_nmf_programs_25.csv")
##jaccard similarity of all NMF programs##
dir<-dir()
P<-dir[grep("nmf_programs_25.csv",dir)]
NMF<-list()
for (i in 1:length(P))
{NMF[[i]]<-read.csv(file=P[[i]],header=T,row.names=1)}
NMF<-do.call("cbind",NMF)
inter_intersect <- apply(NMF , 2, function(x) apply(NMF , 2, function(y) length(intersect(x,y))))
inter_union <- apply(NMF , 2, function(x) apply(NMF , 2, function(y) length(union(x,y))))
inter_intersect<-reshape2::melt(inter_intersect)
inter_union<-reshape2::melt(inter_union)
inter_intersect$union<-inter_union$value
inter_intersect$jaccard<-inter_intersect$value/inter_intersect$union
intersect<-matrix(inter_intersect$jaccard,nrow=ncol(NMF),ncol=ncol(NMF),byrow=F)
colnames(intersect)<-unique(inter_intersect$Var1)
rownames(intersect)<-unique(inter_intersect$Var1)
##NMF module##
BM<-readRDS("BM_MGI_RNA.rds")
data<-BM[unique(NMF$gene),]@assays$RNA@data
data<-apply(data,1,scale)
cor<-cor(data,method="pearson")
dist<-as.dist(1-cor)
p=hclust(dist,method="average")
gene<-as.data.frame(p$label)
gene<-gene[p$order,]
write.csv(gene,"BM_BMt__total_NMF_signature_module_hclust_total_25.csv")
pdf("BM_BMt_ccancertype_total_nmf_module_corrplot_dist_total_25.pdf",width=12,height=5)
dist<-as.dist(1-cor)
plot(p,hang=-1)
dev.off()
