setwd("~/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Real Data/AIRE Deletion/")
load("GSE188870pv.RData")
load("GSE222285pv.RData")
load("GSE224247pv.RData")
load("GSE151012pv.RData")
load("GSE188870mTEChit_mat.RData")
load("genenames.Rdata")

gene_names <- intersect(names(GSE222285mTEChipv),names(GSE224247mTEChipv))
gene_names <- intersect(gene_names, names(GSE151012mTEChipv))
gene_names <- intersect(gene_names, names(GSE188870mTEChipv))
gene_names <- intersect(gene_names, gene.names)

p.mat <- matrix(data = NA, nrow = length(gene.names), ncol = 3)
rownames(p.mat) <- gene.names

p.mat[rownames(p.mat),1] <- GSE222285mTEChipv[rownames(p.mat)]
p.mat[rownames(p.mat),2] <- GSE224247mTEChipv[rownames(p.mat)]
p.mat[rownames(p.mat),3] <- GSE151012mTEChipv[rownames(p.mat)]

X <- GSE188870mTEChit_mat[rownames(p.mat),]
X_spline <- do.call(cbind, lapply(1:ncol(X), function(i) splines::bs(X[, i], df = 6)) )
X_list <- rep(list(   X_spline    ),3)

save(p.mat, X_list, X, X_spline,
     file = "~/Library/CloudStorage/GoogleDrive-ninhtran021998@gmail.com/My Drive/PhD/Replicability/Paper Simulations/Covariate-Assisted/Second Round Codes/RealData2/AIREKO_data.Rdata")
