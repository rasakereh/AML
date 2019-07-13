require(GEOquery)
require(dplyr)
require(udpipe)
require(limma)

setwd("~/Desktop/sharif/9702/Bio/project")

ACCESSION_NO <- "GSE48558"
PLATFORM <- "GPL6244"

gset <- getGEO(ACCESSION_NO, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir="data/")
if (length(gset) > 1) idx <- grep(PLATFORM, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

GENE.SYMBOLS <- data.frame(Gene.symbol=fData(gset)["Gene symbol"])
GENE.SYMBOLS <- GENE.SYMBOLS %>% mutate(probe_id=rownames(GENE.SYMBOLS))
GENE.TITLES <- data.frame(Gene.title=fData(gset)["Gene title"])
GENE.TITLES <- GENE.TITLES %>% mutate(probe_id=rownames(GENE.TITLES))

categories <- as.factor(gset$`phenotype:ch1`)
subpopulations <- as.factor(gset$source_name_ch1)

expressionData <- data.frame(exprs(gset) %>% na.omit())
T.expr <- expressionData %>%
  select(expressionData[which(grepl("T",subpopulations))] %>% colnames())
T.cats <- categories[which(grepl("T",subpopulations))]
B.expr <- expressionData %>%
  select(expressionData[which(grepl("B",subpopulations))] %>% colnames())
B.cats <- categories[which(grepl("B",subpopulations))]

tTGenerator <- function(fl, exprData)
{
  design <- model.matrix(~ fl + 0, exprData)
  colnames(design) <- levels(fl)
  fit <- lmFit(exprData, design)
  cont.matrix <- makeContrasts(Normal-Leukemia, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=500)
  if(!sum(grepl("Gene.symbol", colnames(tT))))
  {
    tT <- tT %>% mutate(probe_id=rownames(tT)) %>% left_join(GENE.SYMBOLS, by="probe_id") %>%
      left_join(GENE.TITLES, by="probe_id")
    rownames(tT) <- tT$probe_id
  }
  tT <- subset(tT, select=c("logFC", "adj.P.Val","P.Value","Gene.symbol","Gene.title"))
  
  return(tT)
}

tT.total <- tTGenerator(categories, expressionData)
tT.B <- tTGenerator(B.cats, B.expr)
tT.T <- tTGenerator(T.cats, T.expr)

total.ups <- (tT.total %>% filter(adj.P.Val < 1e-2 & logFC > 1))[['Gene.symbol']] %>%
  strsplit("///") %>% unlist() %>% unique()
total.downs <- (tT.total %>% filter(adj.P.Val < 1e-2 & logFC < -1))[['Gene.symbol']] %>%
  strsplit("///") %>% unlist() %>% unique()
B.ups <- (tT.B %>% filter(adj.P.Val < 1e-2 & logFC > 1))[['Gene.symbol']] %>%
  strsplit("///") %>% unlist() %>% unique()
B.downs <- (tT.B %>% filter(adj.P.Val < 1e-2 & logFC < -1))[['Gene.symbol']] %>%
  strsplit("///") %>% unlist() %>% unique()
T.ups <- (tT.T %>% filter(adj.P.Val < 1e-2 & logFC > 1))[['Gene.symbol']] %>%
  strsplit("///") %>% unlist() %>% unique()
T.downs <- (tT.T %>% filter(adj.P.Val < 1e-2 & logFC < -1))[['Gene.symbol']] %>%
  strsplit("///") %>% unlist() %>% unique()

write(total.ups, file="results/genes/TotalUps.txt")
write(total.downs, file="results/genes/TotalDowns.txt")
write(B.ups, file="results/genes/BUps.txt")
write(B.downs, file="results/genes/BDowns.txt")
write(T.ups, file="results/genes/TUps.txt")
write(T.downs, file="results/genes/TDowns.txt")

