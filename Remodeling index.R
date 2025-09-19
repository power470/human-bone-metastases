BM<-readRDS("BM_BMT_meta.rds")
scRNA<-readRDS("BM_P_meta.rds")
meta_BMT<-BM@meta.data[c("cancertype","site","celltype")]
meta_P<-scRNA@meta.data[c("cancertype","site","celltype")]
metadat<-rbind(meta_BMT,meta_P)
metadat$RItype[ metadat$celltype%in%c("T_cell","B_cell","Myeloid_cell","Plasma_cell","Osteoclast") ] <- "Immune"
metadat$RItype[ metadat$celltype%in%c("Endothelail_cell","Stromal_cell") ] <- "Stromal"
metadat$RItype[ metadat$celltype%in%c("Tumor_cell") ] <- "Tumor"
# Step2 : Proportion
sample<-c("LUAD","LUSC","LIHC","CHOL","KIRC","BRCA","PRAD","PAAD","COAD","STAD","NPC")
dat01 <- metadat %>% filter(RItype!="Other", site=="BMt")
dat01 <- as.data.frame.matrix(table(dat01$RItype, dat01$cancertype))%>% dplyr::select(sample)
dat01.prop <- t(t(dat01)/colSums(dat01))*100
colnames(dat01.prop) <- paste0("BMT_",colnames(dat01.prop))

dat02 <- metadat %>% filter(RItype!="Other", site=="primary")
dat02 <- as.data.frame.matrix(table(dat02$RItype, dat02$cancertype)) %>% dplyr::select(sample)
dat02.prop <- t(t(dat02)/colSums(dat02))*100
colnames(dat02.prop) <- paste0("PT_",colnames(dat02.prop))

# Step3 : Calculate RI
identical(rownames(dat01),rownames(dat02)) # TRUE
dat <- as.data.frame( t(cbind(dat01.prop, dat02.prop)) )

distval <- philentropy::distance(dat, method = "soergel")
rownames(distval) = colnames(distval) <- rownames(dat)
distval
