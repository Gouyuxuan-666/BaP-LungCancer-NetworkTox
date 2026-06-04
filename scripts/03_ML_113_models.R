#install.packages(c("seqinr", "plyr", "openxlsx", "randomForestSRC", "glmnet", "RColorBrewer"))
#install.packages(c("ade4", "plsRcox", "superpc", "gbm", "plsRglm", "BART", "snowfall"))
#install.packages(c("caret", "mboost", "e1071", "BART", "MASS", "pROC", "xgboost"))

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("mixOmics")
#BiocManager::install("survcomp")
#BiocManager::install("ComplexHeatmap")


#���ð�
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)


#���ù���Ŀ¼
setwd("C:\\Users\\1\\Desktop\\11")
source("refer.ML(1).R")

#��ȡѵ����������ļ�
Train_data <- read.table("data.train.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr=Train_data[,1:(ncol(Train_data)-1),drop=F]
Train_class=Train_data[,ncol(Train_data),drop=F]

#��ȡ������������ļ�
Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors = F)
Test_expr=Test_data[,1:(ncol(Test_data)-1),drop=F]
Test_class=Test_data[,ncol(Test_data),drop=F]
Test_class$Cohort=gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
Test_class=Test_class[,c("Cohort", "Type")]

#��ȡѵ����Ͳ�����Ľ�������
comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[,comgene])
Test_expr <- as.matrix(Test_expr[,comgene])
Train_set = scaleData(data=Train_expr, centerFlags=T, scaleFlags=T)
names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort))
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)

# Clean gene names: strip " /// ..." suffixes from ambiguous Affymetrix probes
colnames(Train_set) <- gsub(" /// .*", "", colnames(Train_set))
colnames(Test_set) <- gsub(" /// .*", "", colnames(Test_set))
colnames(Train_set) <- make.names(colnames(Train_set), unique = TRUE)
colnames(Test_set) <- make.names(colnames(Test_set), unique = TRUE)

# Pre-filter: Lasso screening to reduce 1000+ genes to manageable set
# Stepglm[both] on 1000 variables is impractically slow (exponential complexity)
if(ncol(Train_set) > 200) {
  set.seed(2025)
  cvfit <- cv.glmnet(x = Train_set, y = Train_class[["Type"]], family = "binomial", alpha = 1, nfolds = 5)
  lasso_coef <- coef(cvfit, s = "lambda.1se")
  lasso_genes <- rownames(lasso_coef)[which(lasso_coef[,1] != 0)]
  lasso_genes <- setdiff(lasso_genes, "(Intercept)")
  if(length(lasso_genes) < 10) {
    lasso_coef <- coef(cvfit, s = "lambda.min")
    lasso_genes <- rownames(lasso_coef)[which(lasso_coef[,1] != 0)]
    lasso_genes <- setdiff(lasso_genes, "(Intercept)")
  }
  if(length(lasso_genes) < 5) lasso_genes <- colnames(Train_set)[1:min(100, ncol(Train_set))]
  Train_set <- Train_set[, lasso_genes, drop = FALSE]
  Test_set <- Test_set[, lasso_genes, drop = FALSE]
  cat(sprintf("Pre-filter by Lasso: %d genes retained for ML modeling
", length(lasso_genes)))
}

#��ȡ����ѧϰ�������ļ�
methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)
methods=methodRT$Model
methods <- gsub("-| ", "", methods)


#׼������ѧϰģ�͵Ĳ���
classVar = "Type"         #���÷���ı�����
min.selected.var = 5      #������Ŀ����ֵ
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))


######################����ѵ�������ݹ�������ѧϰģ��######################
#����ģ����ϵ�һ�ֻ���ѧϰ����ɸѡ����
preTrain.var <- list()       #���ڱ�����㷨ɸѡ�ı���
set.seed(seed = 123)         #��������
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method,          #����ѧϰ����
                                 Train_set = Train_set,         #ѵ����Ļ����������
                                 Train_label = Train_class,    #ѵ����ķ�������
                                 mode = "Variable",              #ѡ������ģʽ(ɸѡ����)
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)

#����ģ����ϵڶ��ֻ���ѧϰ��������ģ��
model <- list()            #��ʼ��ģ�ͽ���б�
set.seed(seed = 123)       #��������
Train_set_bk = Train_set
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)
  Variable = preTrain.var[[method[1]]]
  Train_set = Train_set_bk[, Variable]
  Train_label = Train_class
  model[[method_name]] <- RunML(method = method[2],       #����ѧϰ����
                                Train_set = Train_set,         #ѵ����ı�������
                                Train_label = Train_label,    #ѵ����ķ�������
                                mode = "Model",                 #ѡ������ģʽ(����ģ��)
                                classVar = classVar)
  
  #���ĳ�ֻ���ѧϰ����ɸѡ���ı���С����ֵ����÷������Ϊ��
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk; rm(Train_set_bk)
#�������л���ѧϰģ�͵Ľ��
saveRDS(model, "model.MLmodel.rds")

#����������߼��ع�ģ��
FinalModel <- c("panML", "multiLogistic")[2]
if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){    #�����߼��ع�ģ�ͼ���ÿ�������ķ������
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial", 
               data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit)
    return(tmp)
  })
}
#���������Զ�����߼��ع�ģ��
saveRDS(logisticmodel, "model.logisticmodel.rds")


#���ݻ������������ÿ�������ķ���÷�
model <- readRDS("model.MLmodel.rds")            #ʹ�ø�������ѧϰģ�͵�������Ϻ�������÷�
#model <- readRDS("model.logisticmodel.rds")     #ʹ�ö�����߼��ع�ģ�ͼ���÷�
methodsValid <- names(model)                     #��������������Ŀ��ȡ��Ч��ģ��
#���ݻ��������Ԥ�������ķ��յ÷�
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], new_data = rbind.data.frame(Train_set,Test_set))
}
riskTab=as.data.frame(t(do.call(rbind, RS_list)))
riskTab=cbind(id=row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)

#���ݻ��������Ԥ����Ʒ�ķ���
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = rbind.data.frame(Train_set,Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) # ��Ҫ�ϲ����Լ�������������Ϣ�ļ������д���
classTab=cbind(id=row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)

#��ȡÿ�ֻ���ѧϰ����ɸѡ���ı���(ģ�ͻ���)
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file="model.genes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#����ÿ��ģ�͵�AUCֵ
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],      #����ѧϰģ��
                                Test_set = Test_set,        #������ı�������
                                Test_label = Test_class,    #������ķ�������
                                Train_set = Train_set,      #ѵ����ı�������
                                Train_label = Train_class,  #ѵ����ķ�������
                                Train_name = "Train",        #ѵ����ı�ǩ
                                cohortVar = "Cohort",        #GEO��id
                                classVar = classVar)         #�������
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab=cbind(Method=row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)


##############################����AUC��ͼ##############################
#׼��ͼ�ε�����
AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)

#����AUC�ľ�ֵ�Ի���ѧϰģ�ͽ�������
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = T)
AUC_mat <- AUC_mat[names(avg_AUC),]
#��ȡ����ģ��(ѵ����+�������AUC��ֵ���)
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

#������ͼע�͵���ɫ
CohortCol <- c("red", "blue")
if(ncol(AUC_mat)>2){
	CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")}
names(CohortCol) <- colnames(AUC_mat)

#����ͼ��
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = AUC_mat,    #AUCֵ�ľ���
                    avg_Cindex = avg_AUC,        #AUC��ֵ
                    CohortCol = CohortCol,       #���ݼ�����ɫ
                    barCol = "steelblue",        #�Ҳ���״ͼ����ɫ
                    cellwidth = cellwidth, cellheight = cellheight,    #��ͼÿ�����ӵĿ��Ⱥ͸߶�
                    cluster_columns = F, cluster_rows = F)      #�Ƿ�����ݽ��о���

#�����ͼ
pdf(file="model.AUCheatmap.pdf", width=cellwidth * ncol(AUC_mat) + 6, height=cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()
