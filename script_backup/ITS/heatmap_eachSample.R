#! /usr/bin/Rscript
# Rscript heatmap.r output/report/comparsion.txt ./
library(devtools)
library(Biobase)
library(dendextend)
library(gplots)
library(preprocessCore)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
fn<-args[1]
dir<-args[2]
ec<-read.csv(fn, header = TRUE, row.names = 1, sep = ",", check.names=F)
ec <- ec[apply(ec[, -1], 1, function(x) !all(x<15)),]	# one of the cell value within row needs to be greater than 6
zscore <- t(apply(ec, 1, function(x) {(x-mean(x))/sd(x)}))

output = paste0(dir, "/OTUs_heatmap.pdf");
pdf(file = output, width = 25, height = 25)
Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250), rep("brown",323))
heatmap.2(as.matrix(zscore[complete.cases(zscore),]), dendrogram="both", col=bluered(256), scale="none", key=T, keysize=0.5, density.info="none", trace="none",cexRow=2 , cexCol=2, RowSideColors=Label, margins = c(15,15), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,15.0), lwid=c(1.5,0.2,2.5,2.5), srtCol = 45)

dev.off()
