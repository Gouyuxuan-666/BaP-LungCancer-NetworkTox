# ============================================
# Full CIBERSORT Pipeline: Deconvolution → Visualization
# Data:    GSE10072 + GSE19804 merged
# Output:  all in alpha/
# ============================================

# -------------------------------
# STEP 1: CIBERSORT Deconvolution
# -------------------------------
library(IOBR)

setwd("../data")  # Adjust to point to the data/ directory

# Copy expression matrix to alpha
file.copy("C:/Users/1/Desktop/6/merge.normalize.txt",
          "C:/Users/1/Desktop/alpha/merge.normalize.txt", overwrite = TRUE)

# Load
cat("Loading expression matrix...\n")
expr <- read.table("merge.normalize.txt", header = TRUE, sep = "\t",
                   row.names = 1, check.names = FALSE)
cat(sprintf("Loaded: %d genes x %d samples\n", nrow(expr), ncol(expr)))

# Anti-log
max_val <- max(expr, na.rm = TRUE)
cat(sprintf("Max value: %.1f → %s\n", max_val,
            if (max_val < 50) "log2 detected, anti-logging" else "linear, skip"))
if (max_val < 50) expr <- 2^expr

# Groups
con  <- grep("_Control", colnames(expr), value = TRUE)
treat <- grep("_Treat", colnames(expr), value = TRUE)
cat(sprintf("Control: %d, Treat: %d\n", length(con), length(treat)))

# Run CIBERSORT
cat("Running CIBERSORT (perm=1000)...\n")
result <- deconvo_tme(eset = expr, method = "cibersort", arrays = TRUE, perm = 1000)
cat(sprintf("Done: %d samples x %d columns\n", nrow(result), ncol(result)))

# Filter P<0.05
pval_col <- grep("P.value|P-value", colnames(result), value = TRUE)[1]
before <- nrow(result)
result <- result[result[[pval_col]] < 0.05, ]
cat(sprintf("P<0.05 filter: %d → %d samples\n", before, nrow(result)))

# Clean
meta_cols <- c("ID", grep("^P[.-]value|^Correlation|^RMSE", colnames(result), value = TRUE))
cell_cols <- setdiff(colnames(result), meta_cols)
cibersort <- result[, cell_cols, drop = FALSE]
rownames(cibersort) <- result$ID
colnames(cibersort) <- gsub("_CIBERSORT|_ABS$", "", colnames(cibersort))
colnames(cibersort) <- gsub("_", " ", colnames(cibersort))

write.table(cibersort, file = "CIBERSORT-Results.txt", sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf("Saved CIBERSORT-Results.txt: %d x %d\n", nrow(cibersort), ncol(cibersort)))

# -------------------------------
# STEP 2: Visualization
# -------------------------------
library(reshape2)
library(ggpubr)
library(corrplot)

rt <- read.table("CIBERSORT-Results.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
con_idx  <- grep("_Control", rownames(rt))
treat_idx <- grep("_Treat", rownames(rt))
conData   <- rt[con_idx, ]
treatData <- rt[treat_idx, ]
conNum    <- nrow(conData)
treatNum  <- nrow(treatData)
data <- t(rbind(conData, treatData))

cat(sprintf("\nVisualization: Control=%d, Treat=%d\n", conNum, treatNum))

# 1. Barplot
pdf("barplot.pdf", width = 14, height = 7)
col <- rainbow(nrow(data), s = 0.7, v = 0.7)
par(las = 1, mar = c(8, 5, 4, 16), mgp = c(3, 0.1, 0), cex.axis = 1.2)
a1 <- barplot(data, col = col, xaxt = "n", yaxt = "n",
              ylab = "Relative Percent", cex.lab = 1.5, border = NA)
a2 <- axis(2, tick = FALSE, labels = FALSE)
axis(2, a2, paste0(a2 * 100, "%"))
par(srt = 0, xpd = TRUE)
rect(xleft = a1[1] - 0.5, ybottom = -0.01, xright = a1[conNum] + 0.5,
     ytop = -0.08, col = "#6699FFFF")
text(a1[conNum] / 2, -0.045, "Control", cex = 1.8)
rect(xleft = a1[conNum] + 0.5, ybottom = -0.01, xright = a1[length(a1)] + 0.5,
     ytop = -0.08, col = "#E6550DFF")
text((a1[length(a1)] + a1[conNum]) / 2, -0.045, "Treat", cex = 1.8)
legend(par('usr')[2] * 0.98, par('usr')[4], legend = rownames(data),
       col = col, pch = 15, bty = "n", cex = 0.85)
dev.off()
cat("Saved: barplot.pdf\n")

# 2. Boxplot
Type <- gsub("(.*)\\_(.*)", "\\2", colnames(data))
data2 <- cbind(as.data.frame(t(data)), Type)
data2 <- melt(data2, id.vars = "Type")
colnames(data2) <- c("Type", "Immune", "Expression")
group <- levels(factor(data2$Type))
bioCol <- c("#6699FFFF", "#E6550DFF", "#0066FF", "#FF0000",
            "#6E568C", "#7CC767", "#223D6C", "#D20A13",
            "#FFD121", "#088247", "#11AA4D")[1:length(group)]

bp <- ggboxplot(data2, x = "Immune", y = "Expression", fill = "Type",
                xlab = "", ylab = "Fraction", legend.title = "Type",
                width = 0.8, palette = bioCol) +
    rotate_x_text(50) +
    stat_compare_means(aes(group = Type),
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "")),
                       label = "p.signif")
pdf("immune.diff.pdf", width = 9, height = 6)
print(bp)
dev.off()
cat("Saved: immune.diff.pdf\n")

# 3. Correlation Heatmap
treatData <- treatData[, apply(treatData, 2, sd) > 0]
pdf("corHeatmap.pdf", width = 12, height = 12)
corrplot(corr = cor(treatData, method = "spearman"),
         method = "color", order = "hclust",
         tl.col = "black", number.cex = 0.8, addCoef.col = "black",
         col = colorRampPalette(c("#6699FFFF", "white", "#E6550DFF"))(50))
dev.off()
cat("Saved: corHeatmap.pdf\n")

# Summary
cat("\n========================================\n")
cat("alpha/ pipeline complete:\n")
cat("  CIBERSORT-Results.txt\n")
cat("  barplot.pdf\n")
cat("  immune.diff.pdf\n")
cat("  corHeatmap.pdf\n")
cat("========================================\n")
