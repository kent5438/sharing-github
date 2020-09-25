library(devtools)
library(Biobase)
library(dendextend)
library(gplots)
library(preprocessCore)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
fn<-args[1]
dir<-args[2]

ec<-read.table(fn, header = TRUE, sep = "\t")
ec$Gene <- NULL

# genMat_norm is median-by-ratio normalization function from Anders et. al., 2010.
# MedianNorm -- Geometric mean (?)
# geomeans <- exp(rowMeans(log(Data)))
# out <- apply(Data, 2, function(cnts) median((cnts/geomeans)[geomeans > 0]))
# GetNormalizedMat
# Out = Data/Sizes (the Sizes is the output of MedianNorm)

row_sub <- apply(ec, 1, function(row) all(row!=0))
ec.no0 = ec[row_sub,]
ec.no0 <- log2(ec.no0)
ec.no0.2 <- normalize.quantiles(as.matrix(ec.no0)) 
rownames(ec.no0.2) <- rownames(ec.no0)
colnames(ec.no0.2) <- colnames(ec.no0)

# HierarchicalClustering
output = paste0(dir, "/HC.png");
dist1 <- dist(t(ec.no0.2))
png(filename = output, width = 1600, height = 1200)
par(col="black")
hclust1 = hclust(dist1)
plot(hclust1, hang=-1)
dev.off()





