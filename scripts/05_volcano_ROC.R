# Step 24-25: Volcano (gene expression) + ROC curves for ML-selected genes
options(repos = "https://cloud.r-project.org")

library(pROC)
library(ggplot2)
library(reshape2)

setwd("C:\\Users\\1\\Desktop\\11")

# Load data
Train_data <- read.table("data.train.txt", header=T, sep="\t", check.names=F, row.names=1)
Train_class <- Train_data[, ncol(Train_data)]
Train_expr <- Train_data[, -ncol(Train_data), drop=FALSE]

Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1)
Test_class <- Test_data[, ncol(Test_data)]
Test_expr <- Test_data[, -ncol(Test_data), drop=FALSE]

# Clean gene names
colnames(Train_expr) <- gsub(" /// .*", "", colnames(Train_expr))
colnames(Train_expr) <- make.names(colnames(Train_expr), unique=TRUE)
colnames(Test_expr) <- gsub(" /// .*", "", colnames(Test_expr))
colnames(Test_expr) <- make.names(colnames(Test_expr), unique=TRUE)

# Load SHAP results
shap_imp <- read.table("23.SHAP.importance.txt", header=T, sep="\t", check.names=F)
shap_genes <- intersect(shap_imp$Feature, colnames(Train_expr))
cat("SHAP genes available:", length(shap_genes), "\n")

# Compute fold change and p-value for SHAP genes using train data
train_ctrl <- Train_expr[Train_class == 0, shap_genes, drop=FALSE]
train_treat <- Train_expr[Train_class == 1, shap_genes, drop=FALSE]

deg_stats <- data.frame(
  Gene = shap_genes,
  meanCtrl = colMeans(train_ctrl),
  meanTreat = colMeans(train_treat),
  logFC = colMeans(train_treat) - colMeans(train_ctrl),
  pvalue = sapply(shap_genes, function(g) {
    tryCatch(t.test(Train_expr[Train_class==1, g], Train_expr[Train_class==0, g])$p.value,
             error=function(e) NA)
  }),
  stringsAsFactors=FALSE
)
deg_stats$negLogP <- -log10(deg_stats$pvalue)
deg_stats$significant <- deg_stats$pvalue < 0.05
deg_stats <- deg_stats[order(deg_stats$negLogP, decreasing=TRUE),]

# Volcano-like plot for SHAP genes
pdf("24.volcano.SHAPgenes.pdf", width=10, height=8)
ggplot(deg_stats, aes(x=logFC, y=negLogP)) +
  geom_point(aes(color=significant), size=4, alpha=0.8) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="grey60")) +
  geom_text(aes(label=Gene), size=3, vjust=-1, hjust=0.5, check_overlap=TRUE) +
  geom_vline(xintercept=0, linetype="dashed", alpha=0.3) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.3) +
  theme_minimal() +
  labs(title="Gene Expression: Treat vs Control (Training Set)",
       subtitle=paste(nrow(deg_stats), "SHAP-selected genes from Lasso+glmBoost"),
       x="Mean Expression Difference (Treat - Control)",
       y="-log10(P-value)") +
  theme(legend.position="none")
dev.off()

# Expression boxplots for top 10 SHAP genes
top_plot_genes <- head(deg_stats$Gene, 10)
plot_data <- data.frame()
for (g in top_plot_genes) {
  plot_data <- rbind(plot_data,
    data.frame(Gene=g, Expression=Train_expr[Train_class==0, g], Group="Control"),
    data.frame(Gene=g, Expression=Train_expr[Train_class==1, g], Group="Treat"))
}

pdf("24.boxplot.SHAPgenes.pdf", width=14, height=10)
ggplot(plot_data, aes(x=Group, y=Expression, fill=Group)) +
  geom_boxplot(outlier.size=1) +
  facet_wrap(~Gene, scales="free_y", ncol=5) +
  scale_fill_manual(values=c("Control"="#4575B4", "Treat"="#D73027")) +
  theme_minimal() +
  labs(title="Expression of Top SHAP Genes: Control vs Treat", y="Expression Level") +
  theme(strip.text=element_text(face="bold"))
dev.off()

# ROC curves for top 5 SHAP genes (single-gene models on test set)
top_genes <- head(deg_stats$Gene, 5)
cat("\nTop 5 genes for ROC:", paste(top_genes, collapse=", "), "\n")

pdf("25.ROC.SHAPgenes.pdf", width=8, height=8)
plot(0:1, 0:1, type="n", xlab="1 - Specificity", ylab="Sensitivity",
     main="ROC Curves: Top SHAP Genes (Test Set: GSE19804)")
abline(a=0, b=1, lty=2, col="grey60")
cols <- c("#D73027", "#4575B4", "#1A9850", "#FDAE61", "#762A83")
aucs <- c()

for (i in seq_along(top_genes)) {
  g <- top_genes[i]
  x_train <- Train_expr[[g]]
  x_test <- Test_expr[[g]]
  fit <- glm(Train_class ~ x_train, family="binomial")
  pred <- predict(fit, newdata=data.frame(x_train=x_test), type="response")

  if (length(unique(Test_class)) >= 2) {
    roc_obj <- roc(Test_class, pred, quiet=TRUE)
    lines(1-roc_obj$specificities, roc_obj$sensitivities, col=cols[i], lwd=2.5)
    aucs <- c(aucs, auc(roc_obj))
    cat(sprintf("%s: Test AUC = %.4f\n", g, auc(roc_obj)))
  }
}
legend("bottomright",
       legend=paste0(top_genes, " (AUC=", round(aucs, 3), ")"),
       col=cols, lwd=2.5, cex=0.9, bty="n")
dev.off()

# Combined multi-gene logistic model
cat("\n--- Combined model ---\n")
common <- intersect(top_genes, colnames(Train_expr))
if (length(common) >= 3) {
  train_df <- as.data.frame(Train_expr[, common, drop=FALSE])
  test_df <- as.data.frame(Test_expr[, common, drop=FALSE])
  train_df$y <- Train_class
  fit_multi <- glm(y ~ ., family="binomial", data=train_df)
  pred_multi <- predict(fit_multi, newdata=test_df, type="response")
  multi_roc <- roc(Test_class, pred_multi, quiet=TRUE)
  cat(sprintf("Multi-gene (5 SHAP) AUC: %.4f\n", auc(multi_roc)))
}

# Save stats
write.table(deg_stats, "24.SHAPgenes.stats.txt", sep="\t", quote=F, row.names=F)

cat("\nVolcano + Boxplot + ROC completed.\n")
