setwd("/Users/ninht/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Real Data/AIRE Deletion/")
library(readxl)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)


x <- read.delim("GSE151012_Deseq_all_results_AireKO_B6.txt")
head(x)
(x$Aire_WT_Aire_KO.pvalue)[which(x$X == "Aire")]
(x$Aire_WT_Aire_KO.Direction)[which(x$X == "Aire")]
x$Aire_WT_Aire_KO.log2FoldChange[which(x$X == "Aire")]
x$Aire_WT_Aire_KO.baseMean[which(x$X == "Aire")]

x <- DGEList(counts = x[,2:14], genes = x$X)
colnames(x$counts)
group <- c("AIRE_het","AIRE_het","AIRE_het","AIRE_het","AIRE_het",
            "AIRE_wt","AIRE_wt","AIRE_wt","AIRE_wt",
            "AIRE_ko","AIRE_ko","AIRE_ko","AIRE_ko")
length(group)
group <- as.factor(group)
x$samples$group <- group
keep.exprs <- filterByExpr(x, group=group, min.total.count = 100)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
contr.matrix <- makeContrasts(
  wtvsko  = AIRE_wt-AIRE_ko,
  levels = colnames(design))
contr.matrix
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

p <- 2*(1-pt(q = abs(efit$t), df = efit$df.total))
p <- c(p)
names(p) <- x$genes$genes
p["Aire"]
GSE151012mTEChipv <- p
GSE151012mTEChitv <- efit$t
GSE151012mTEChizv <- qnorm(pt(q = efit$t, df = efit$df.total))
names(GSE151012mTEChitv) <- names(GSE151012mTEChipv)
names(GSE151012mTEChizv) <- names(GSE151012mTEChipv)
sum(p.adjust(p = GSE151012mTEChipv, method = "BH") <= 0.01)
save(GSE151012mTEChipv, file = "GSE151012pv.RData")
save(GSE151012mTEChitv, file = "GSE151012tv.RData")
save(GSE151012mTEChizv, file = "GSE151012zv.RData")

