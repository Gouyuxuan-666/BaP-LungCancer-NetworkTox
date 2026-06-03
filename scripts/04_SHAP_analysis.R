# Step 23: SHAP analysis for best ML model
library(glmnet)
library(plyr)
library(MASS)
library(mboost)
library(ggplot2)
library(pROC)

source("refer.ML(1).R")
setwd("C:\\Users\\1\\Desktop\\11")

# Load data
Train_data <- read.table("data.train.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr <- as.matrix(Train_data[,1:(ncol(Train_data)-1),drop=F])
Train_class <- Train_data[, ncol(Train_data), drop=F]

Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)
Test_expr <- as.matrix(Test_data[,1:(ncol(Test_data)-1),drop=F])
Test_class <- Test_data[, ncol(Test_data), drop=F]

comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[,comgene])
Test_expr <- as.matrix(Test_expr[,comgene])
Train_set <- scaleData(data=Train_expr, centerFlags=T, scaleFlags=T)
Test_class$Cohort <- gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
Test_class <- Test_class[,c("Cohort", "Type")]
Test_set <- scaleData(data=Test_expr, cohort=Test_class$Cohort, centerFlags=T, scaleFlags=T)

# Clean gene names
colnames(Train_set) <- gsub(" /// .*", "", colnames(Train_set))
colnames(Test_set) <- gsub(" /// .*", "", colnames(Test_set))
colnames(Train_set) <- make.names(colnames(Train_set), unique=TRUE)
colnames(Test_set) <- make.names(colnames(Test_set), unique=TRUE)

# Load best model
model <- readRDS("model.MLmodel.rds")
methodsValid <- names(model)

# Find best model by AUC
AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1)
avg_AUC <- rowMeans(AUC_mat, na.rm=TRUE)
best_method <- names(which.max(avg_AUC))
cat("Best model:", best_method, "Avg AUC:", round(max(avg_AUC, na.rm=TRUE), 4), "\n")

# Extract features from best model
best_fit <- model[[best_method]]
best_features <- ExtractVar(best_fit)
cat("Selected features (", length(best_features), "):", paste(best_features, collapse=", "), "\n")

# Use combined data for SHAP
X_all <- rbind(Train_set[, best_features, drop=F], Test_set[, best_features, drop=F])
y_all <- c(Train_class[["Type"]], Test_class[["Type"]])

# SHAP using fastshap (model-agnostic)
options(repos = "https://cloud.r-project.org")
if (!requireNamespace("fastshap", quietly=TRUE)) install.packages("fastshap")
library(fastshap)

# Build explainer using logistic regression on best features
explain_data <- as.data.frame(Train_set[, best_features, drop=F])
explain_data$Type <- Train_class[["Type"]]
final_glm <- glm(Type ~ ., family="binomial", data=explain_data)

# Prediction function for SHAP
pfun <- function(object, newdata) {
  predict(object, newdata=as.data.frame(newdata), type="response")
}

# Compute SHAP values
set.seed(123)
X_explain <- as.data.frame(rbind(Train_set[, best_features, drop=F],
                                  Test_set[, best_features, drop=F]))
shap_values <- explain(object=final_glm, X=X_explain, pred_wrapper=pfun, nsim=100)

# SHAP summary bar plot
shap_importance <- colMeans(abs(shap_values))
shap_importance <- sort(shap_importance, decreasing=TRUE)
pdf("23.SHAP.barplot.pdf", width=8, height=6)
par(mar=c(5,10,3,2))
barplot(rev(shap_importance), horiz=TRUE, las=2, col="steelblue",
        main="SHAP Feature Importance", xlab="Mean |SHAP|", border=NA)
dev.off()

# SHAP beeswarm plot
for (pkg in c("ggplot2", "reshape2")) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)
}
library(ggplot2)
library(reshape2)

shap_long <- melt(shap_values)
colnames(shap_long) <- c("Sample", "Feature", "SHAP")
shap_long$Feature <- factor(shap_long$Feature, levels=names(shap_importance))

pdf("23.SHAP.beeswarm.pdf", width=10, height=6)
ggplot(shap_long, aes(x=SHAP, y=Feature)) +
  geom_jitter(aes(color=SHAP), width=0, height=0.2, alpha=0.6, size=1.5) +
  scale_color_gradient2(low="blue", mid="grey90", high="red", midpoint=0) +
  theme_minimal() +
  labs(title=paste("SHAP Beeswarm Plot -", best_method), x="SHAP Value") +
  theme(legend.position="right")
dev.off()

# SHAP dependence plots for top features
top_features <- names(shap_importance)[1:min(5, length(shap_importance))]
for (feat in top_features) {
  pdf(paste0("23.SHAP.dependence.", feat, ".pdf"), width=7, height=5)
  plot(X_explain[[feat]], shap_values[, feat],
       xlab=paste(feat, "Expression"), ylab="SHAP Value",
       main=paste("SHAP Dependence:", feat), pch=19, col=rgb(0,0,0,0.4))
  lines(lowess(X_explain[[feat]], shap_values[, feat]), col="red", lwd=2)
  dev.off()
}

# Save SHAP values
write.table(shap_values, "23.SHAP.values.txt", sep="\t", quote=F, row.names=TRUE, col.names=NA)
write.table(data.frame(Feature=names(shap_importance), MeanAbsSHAP=shap_importance),
            "23.SHAP.importance.txt", sep="\t", quote=F, row.names=FALSE)

cat("SHAP analysis completed.\n")
