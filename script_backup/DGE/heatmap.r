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
ec<-read.table(fn, header = TRUE, sep = "\t")
ec$X <-NULL # There is a extra tab in each line

rownames(ec) <-ec$Gene
ec$Gene <- NULL
output = paste0(dir, "/heatmap.png");
png(filename = output, width = 1600, height = 1200)
Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250), rep("brown",323))
#heatmap.2(log2(as.matrix(ec)), col=redblue(256), dendrogram="both", scale="none", key=T, keysize=0.5, density.info="none", trace="none",cexCol=1.2, labRow=NA, RowSideColors=Label, margins = c(10,10), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0), lwid=c(1.5,0.2,2.5,2.5))
heatmap.2(log2(as.matrix(ec)), col=bluered(256), dendrogram="both", scale="none", key=T, keysize=0.5, density.info="none", trace="none",cexCol=1.2, labRow=NA, RowSideColors=Label, margins = c(10,10), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0), lwid=c(1.5,0.2,2.5,2.5))
dev.off()
