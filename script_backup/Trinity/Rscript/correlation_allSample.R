#! /usr/bin/Rscript

library(corrplot)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
fn<-args[1]
output<-args[2]

ec<-read.table(fn, header = TRUE, row.names = 1, sep = "\t", check.names=F)
ec <- ec[apply(ec[,-1], 1, function(x) !all(x<15)),]    # one of the cell value within row needs to be greater than 15
zscore <- t(apply(ec, 1, function(x) {(x-mean(x))/sd(x)}))
data <- cor(zscore)

png(file = output, width = 10, height = 10, units = "in", res = 300, pointsize = 15)
par(cex=0.7)
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
#corrplot(data, method="color", col=col(200), order="hclust", addCoef.col = "black", tl.col="black", tl.srt=45, sig.level = 0.01, insig = "label_sig", diag=FALSE, tl.cex=1/par("cex"), cl.cex=1/par("cex"))
corrplot(data, method="color", order="hclust", col=col(200), addCoef.col="black", tl.col="black", tl.srt=45, tl.cex=1/par("cex"), cl.cex=1/par("cex"))
dev.off()
