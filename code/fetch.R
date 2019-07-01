require(GEOquery)
require(dplyr)
require(ggplot2)
require(pheatmap)
require(tsne)

setwd("~/Desktop/sharif/9702/Bio/project")

ACCESSION_NO <- "GSE48558"
PLATFORM <- "GPL6244"

gset <- getGEO(ACCESSION_NO, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir="data/")
if (length(gset) > 1) idx <- grep(PLATFORM, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gsms <- paste0("00000000000002222232322322332222233232364236423626",
               "62662662262262262272162242455721621162116721672168",
               "33683368336322222222232222222333333334444444001000",
               "77777771555541111111")
categories <- c("AML Primary", "Normal T", "AML B", "AML T", "Normal G", "Normal Mono",
            "AML Cell", "Normal B", "Normal CD34")

for(i in 1:length(categories))
{
  gsms <- gsub(as.character(i-1), paste0(categories[i], "-"), gsms)
}
groups <- strsplit(gsms, '-')[[1]]

expressionData <- exprs(gset) %>% na.omit()
centeredExpression <- expressionData %>% t() %>% scale(scale=F) %>% t()
resPCA <- prcomp(centeredExpression)
resTSNE <- data.frame(tsne(centeredExpression %>% t(), max_iter=200), Category=groups)
# (D)ata * (R)otation = (X)
# D: 4577*170, rows: g_i, cols: d_i
# R: 170*170, rows: d_i, cols: pc_i
# X: 4577*170, rows: g_i, cols: pc_i

groupedPCA <- data.frame(resPCA$rotation, Category=groups)

PCA.T.compare <- groupedPCA %>% filter(Category=="Normal T" | Category=="AML T")
PCA.B.compare <- groupedPCA %>% filter(Category=="Normal B" | Category=="AML B")
PCA.total.compare <- groupedPCA %>% mutate(row.num=1:dim(groupedPCA)[1]) %>%
                mutate(Category=ifelse(row.num %in% grep("AML",Category), "AML", "Normal")) %>% select(-row.num)
SNE.T.compare <- resTSNE %>% filter(Category=="Normal T" | Category=="AML T")
SNE.B.compare <- resTSNE %>% filter(Category=="Normal B" | Category=="AML B")
SNE.total.compare <- resTSNE %>% mutate(row.num=1:dim(groupedPCA)[1]) %>%
  mutate(Category=ifelse(row.num %in% grep("AML",Category), "AML", "Normal")) %>% select(-row.num)
ggplot(PCA.T.compare, aes(PC1, PC2, color=Category)) + geom_point(size=2.8) + theme_bw()
ggplot(PCA.B.compare, aes(PC1, PC2, color=Category)) + geom_point(size=2.8) + theme_bw()
ggplot(PCA.total.compare, aes(PC1, PC2, color=Category)) + geom_point(size=2.8) + theme_bw()
ggplot(SNE.T.compare, aes(X1, X2, color=Category)) + geom_point(size=2.8) + theme_bw()
ggplot(SNE.B.compare, aes(X1, X2, color=Category)) + geom_point(size=2.8) + theme_bw()
ggplot(SNE.total.compare, aes(X1, X2, color=Category)) + geom_point(size=2.8) + theme_bw()
pdf("results/heatmap", height=30, width=30)
pheatmap(cor(expressionData))
dev.off()
