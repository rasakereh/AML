require(GEOquery)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(plotly)
require(Rtsne)

setwd("~/Desktop/sharif/9702/Bio/project")

ACCESSION_NO <- "GSE48558"
PLATFORM <- "GPL6244"

gset <- getGEO(ACCESSION_NO, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir="data/")
if (length(gset) > 1) idx <- grep(PLATFORM, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

category <- as.factor(gset$`phenotype:ch1`)
subpopulation <- as.factor(gset$source_name_ch1)

expressionData <- exprs(gset) %>% na.omit()
centeredExpression <- expressionData %>% t() %>% scale(scale=F) %>% t()
resPCA <- prcomp(centeredExpression)
resTSNE <- Rtsne(centeredExpression %>% t(), max_iter=300)
# (D)ata * (R)otation = (X)
# D: 4577*170, rows: g_i, cols: d_i
# R: 170*170, rows: d_i, cols: pc_i
# X: 4577*170, rows: g_i, cols: pc_i

groupedPCA <- data.frame(resPCA$rotation, Category=category, Subpop=subpopulation)
groupedTSNE <- data.frame(resTSNE$Y, Category=category, Subpop=subpopulation)

PCA.T.compare <- groupedPCA %>% filter(grepl("T",Subpop))
PCA.B.compare <- groupedPCA %>% filter(grepl("B",Subpop))
SNE.T.compare <- groupedTSNE %>% filter(grepl("T",Subpop))
SNE.B.compare <- groupedTSNE %>% filter(grepl("B",Subpop))

pdf("results/BasicVisualisation.pdf")
ggplot(PCA.T.compare, aes(PC1, PC2, color=Category)) + geom_point(size=2.8) +
  labs(title="T Cells(PCA)") + theme_bw()
ggplot(PCA.B.compare, aes(PC1, PC2, color=Category)) + geom_point(size=2.8) +
  labs(title="B Cells(PCA)") + theme_bw()
ggplot(groupedPCA, aes(PC1, PC2, color=Category)) + geom_point(size=2.8) +
  labs(title="Mixture of Cells(PCA)") + theme_bw()
ggplot(SNE.T.compare, aes(X1, X2, color=Category)) + geom_point(size=2.8) +
  labs(title="T Cells(t-SNE)") + theme_bw()
ggplot(SNE.B.compare, aes(X1, X2, color=Category)) + geom_point(size=2.8) +
  labs(title="B Cells(t-SNE)") + theme_bw()
ggplot(groupedTSNE, aes(X1, X2, color=Category)) + geom_point(size=2.8) +
  labs(title="Mixture of Cells(t-SNE)") + theme_bw()
dev.off()
pdf("results/boxplot.pdf", height=20, width=40)
boxplot(expressionData, names=subpopulation)
dev.off()
pdf("results/heatmap.pdf", height=30, width=30)
pheatmap(cor(expressionData), labels_row=category, labels_col=category)
dev.off()
