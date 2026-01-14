setwd("/Users/ninht/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Real Data/AIRE Deletion/")
library(readxl)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

x <- read.delim("GSE188870.txt")
x <- x[-which(is.na(x$Symbol)),]
x <- DGEList(counts = x[,4:21], genes = x$Symbol)
#genes <- select(x = org.Mm.eg.db, keys=geneid, columns=c("SYMBOL","ENTREZID"),
#                keytype="ENTREZID")
head(genes)
group <- c("mTEChi_wild",
           "mTEClo_wild",
           "cTEC_wild",
           "mTEChi_wild",
           "mTEClo_wild",
           "cTEC_wild",
           "mTEChi_wild",
           "mTEClo_wild",
           "cTEC_wild",
           "mTEChi_KO",
           "mTEClo_KO",
           "cTEC_KO",
           "mTEChi_KO",
           "mTEClo_KO",
           "cTEC_KO",
           "mTEChi_KO",
           "mTEClo_KO",
           "cTEC_KO")
group <- as.factor(group)
x$samples$group <- group

keep.exprs <- filterByExpr(x, group=group, min.total.count = 10)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
contr.matrix <- makeContrasts(
  cTECwildvsKO = cTEC_wild-cTEC_KO,
  mTEChiwildevsK0 = mTEChi_wild-mTEChi_KO,
  mTEClowildevsK0  = mTEClo_wild-mTEClo_KO,
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
two.sided.p <- efit$p.value
rownames(p) <- c(x$genes[[1]])
GSE188870mTEChipv <- p[,2]
GSE188870mTEChitv <- efit$t[,2]
GSE188870mTEChit_mat <- efit$t
rownames(GSE188870mTEChit_mat) <- names(GSE188870mTEChipv)
GSE188870mTEChizv <- qnorm(pt(q = efit$t[,2], df = efit$df.total))
names(GSE188870mTEChitv) <- names(GSE188870mTEChipv)
names(GSE188870mTEChizv) <- names(GSE188870mTEChipv)
save(GSE188870mTEChipv, file = "GSE188870pv.RData")
save(GSE188870mTEChitv, file = "GSE188870tv.RData")
save(GSE188870mTEChizv, file = "GSE188870zv.RData")
save(GSE188870mTEChit_mat, file = "GSE188870mTEChit_mat.RData")

GSE188870mTEChipv["Aire"]
sum(p.adjust(p = GSE188870mTEChipv, method = "BH") <= 0.01)

GSE188870mTEChit <- efit$t[,2]
names(GSE188870mTEChit) <- c(x$genes[[1]])
