#! /usr/bin/Rscript

library(corrplot)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
fn<-args[1]
output<-args[2]

ec<-read.table(fn, header = TRUE, row.names = 1, sep = "\t", check.names=F, quote="\"")
data <- cor(ec)

pdf(file = output, width = 25, height = 25, pointsize = 15)
par(cex=0.7)
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
corrplot(data, method="circle", col=col(200), type="upper", order="hclust", addCoef.col = "black", tl.col="black", tl.srt=45, sig.level = 0.01, insig = "blank", diag=FALSE, tl.cex=1/par("cex"), cl.cex=1/par("cex"))
dev.off()
