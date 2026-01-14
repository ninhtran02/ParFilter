setwd("/Users/ninht/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Real Data/AIRE Deletion/")
library(readxl)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

x <- read.csv("GSE222285.csv")
geneid <- x$X
genes <- select(x = org.Mm.eg.db, keys=geneid, columns=c("SYMBOL","ENTREZID"),
                keytype="ENSEMBL")
genes <- genes[!duplicated(genes$ENSEMBL),]
x <- DGEList(counts = x[-which(is.na(genes$SYMBOL)),2:24], genes = genes[-which(is.na(genes$SYMBOL)),])

group <- c("mTEClo_WT", "mTEClo_WT", "mTEClo_WT", "mTEClo_WT", "mTEClo_WT", "mTEClo_WT", "mTEClo_WT", "mTEClo_WT",
           "mTEClo_AireKO", "mTEClo_AireKO", "mTEClo_AireKO", "mTEClo_AireKO",
           "mTEChi_WT", "mTEChi_WT", "mTEChi_WT", "mTEChi_WT", "mTEChi_WT", "mTEChi_WT", "mTEChi_WT",
           "mTEChi_AireKO", "mTEChi_AireKO", "mTEChi_AireKO", "mTEChi_AireKO")
length(group)
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
  WT_vs_KO = mTEChi_WT - mTEChi_AireKO,
  levels = colnames(design))
contr.matrix
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
dim(efit$p.value)

p <- 2*(1-pt(q = abs(efit$t), df = efit$df.total))
names(p) <- c(x$genes$SYMBOL)
p["Aire"]
GSE222285mTEChipv <- p
GSE222285mTEChitv <- efit$t
GSE222285mTEChizv <- qnorm(pt(q = efit$t, df = efit$df.total))
names(GSE222285mTEChitv) <- names(GSE222285mTEChipv)
names(GSE222285mTEChizv) <- names(GSE222285mTEChipv)
save(GSE222285mTEChipv, file = "GSE222285pv.RData")
save(GSE222285mTEChitv, file = "GSE222285tv.RData")
save(GSE222285mTEChizv, file = "GSE222285zv.RData")
sum(p.adjust(p = GSE222285mTEChipv, method = "BH") <= 0.01)


