setwd("/Users/ninht/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Real Data/AIRE Deletion/")
library(readxl)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

x <- read.delim("GSE224247.txt", header = TRUE)
geneid <- x$X
x <- DGEList(counts = x[,8:13], genes = geneid)
dim(x$counts)
head(x$counts)
group <- c("WT","WT","WT","KO","KO","KO")
length(group)
group <- as.factor(group)
x$samples$group <- group

keep.exprs <- filterByExpr(x, group=group, min.total.count = 10)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
#x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
contr.matrix <- makeContrasts(
  WT_vs_KO = WT - KO,
  levels = colnames(design))
contr.matrix
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
p <- 2*(1-pt(q = abs(efit$t[,1]), df = efit$df.total))
sum(p.adjust(p = p, method = "BH") <= 0.05)
names(p) <- x$genes$genes
GSE224247mTEChipv <- p
GSE224247mTEChitv <- efit$t
GSE224247mTEChizv <- qnorm(pt(q = efit$t, df = efit$df.total))
names(GSE224247mTEChitv) <- names(GSE224247mTEChipv)
names(GSE224247mTEChizv) <- names(GSE224247mTEChipv)
sum(p.adjust(p = GSE224247mTEChipv , method = "BH") <= 0.01)
p["Aire"]
save(GSE224247mTEChipv, file = "GSE224247pv.RData")
save(GSE224247mTEChitv, file = "GSE224247tv.RData")
save(GSE224247mTEChizv, file = "GSE224247zv.RData")
