library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(relaimpo)
cellratio$sex<-as.factor(cellratio$sex)
cellratio$treatent<-as.factor(cellratio$treatent)
cellratio$cancertype<-as.factor(cellratio$cancertype)
dependent_vars <- colnames(cellratio)  # 多个因变量
independent_vars <- c("cancertype","treatent","sex")
cellratio_1<-cellratio[,!(colnames(cellratio)%in%c("cancertype","treatent","site","sample"))] # 多个因变量
independent_vars <- c("cancertype","site","treatent")
calc_relative_importance <- function(dep_var) {
  formula <- as.formula(paste(dep_var, "~", paste(independent_vars, collapse = "+")))
  model <- lm(formula, data = cellratio)
  calc.relimp(model, type = "lmg", rela = TRUE)
}
results <- lapply(dependent_vars, calc_relative_importance)
