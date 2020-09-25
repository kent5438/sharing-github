# Rscript heatmap.r output/report/comparsion.txt ./
library(devtools)
library(Biobase)
library(dendextend)
library(gplots)
library(preprocessCore)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
intersectfn<-args[1]
unionfn<-args[2]
dir<-args[3]

#### intersect 

ec<-read.table(intersectfn, header = TRUE, sep = "\t")
ec$X <-NULL # There is a extra tab in each line

names(ec) <-gsub(".genMat.result", "", names(ec))
rownames(ec) <-ec$Gene
ec$Gene <- NULL
output = paste0(dir, "/intersect-heatmap.pdf");
pdf(file = output, width = 25, height = 25)
Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250), rep("brown",323))
heatmap.2(log2(as.matrix(ec)), col=bluered(256), dendrogram="both", scale="none", key=T, keysize=0.5, density.info="none", trace="none",cexCol=1.5,cexRow=1.5 ,RowSideColors=Label, margins = c(10,10), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,15.0), lwid=c(1.5,0.2,2.5,2.5), srtCol = 45)
dev.off()

#### union
ec<-read.table(unionfn, header = TRUE, row.names = 1, sep = "\t")
ec$X <-NULL # There is a extra tab in each line
names(ec) <-gsub(".genMat.result", "", names(ec))

output = paste0(dir, "/union-heatmap.pdf");
pdf(file = output, width = 25, height = 25)
Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250), rep("brown",323))
heatmap.2(log2(as.matrix(ec)), col=redblue(256), dendrogram="both", scale="none", key=T, keysize=0.5, density.info="none", trace="none",cexCol=1.2, labRow=NA, RowSideColors=Label, margins = c(10,10), lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,15.0), lwid=c(1.5,0.2,2.5,2.5), srtCol = 45)
dev.off()

