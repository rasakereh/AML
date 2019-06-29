require(GEOquery)
require(dplyr)
require(ggplot2)

setwd("~/Desktop/sharif/9702/Bio/project")

ACCESSION_NO <- "GSE48558"
PLATFORM <- "GPL6244"

gset <- getGEO(ACCESSION_NO, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir="data/")
if (length(gset) > 1) idx <- grep(PLATFORM, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gsms <- paste0("00000000000000000000000000000000000000001000100000",
               "00000000000000000010100010111101001100110101010101",
               "00010001000000000000000000000000000001111111001000",
               "11111111111111111111")
groups <- strsplit(gsub("1", "Healthy-", gsub("0", "AML-", gsms)), '-')[[1]]

expressionData <- exprs(gset) %>% na.omit()
centeredExpression <- expressionData %>% t() %>% scale(scale=F) %>% t()
resPCA <- prcomp(centeredExpression)
# (D)ata * (R)otation = (X)
# D: 4577*170, rows: g_i, cols: d_i
# R: 170*170, rows: g_i, rows: pc_i
# X: 4577*170, rows: pc_i, rows: d_i

groupedPCA <- data.frame(resPCA$rotation, Category=groups)

ggplot(groupedPCA, aes(PC1, PC2, color=Category)) + geom_point(size=2.8) + theme_bw()
pdf("results/heatmap", height=30, width=30)
pheatmap(cov(resPCA$rotation %>% t()))
dev.off()
