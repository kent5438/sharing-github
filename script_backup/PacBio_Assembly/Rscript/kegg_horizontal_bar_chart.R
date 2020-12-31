# Rscript /export/EC1680U/Rscript/Pacbio/kegg_horizontal_bar_chart.R 4_Annotation/KEGG/PB17012_NAI5_ec2kegg.txt 4_Annotation/KEGG ead
# Author: kentchen
# Date: 2017-06-12

library(ggplot2)

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
df <- args[1]
dir <- args[2]
code <- args[3]

ylabel = paste0("log2 of EC shared with Ref [", code)
ylabel = paste0(ylabel, "]")

ec <- read.table(df, header=T, sep="\t", check.names=F)
#ec.filter <- ec[(ec$P.value)<0.01,]	# filter pvalue < 0.01
ec <- ec[, c(2, 3, 7)]
colnames(ec) <- c("PathwayName", "Category", "EC_from_ref")
ec$EC_from_ref <- log2(ec$EC_from_ref)
ec <- ec[ec$EC_from_ref>0,]	# discard 0 hit count of the 8th column
ec <- ec[order(ec$Category),]
ec$PathwayName <- factor(ec$PathwayName, as.character(ec$PathwayName))
output = paste0(dir, "/kegg_barchart.png")

png(output, units="in" , width=25, height=20, res=150)
#ggplot(ec, aes(x=PathwayName, y=EC_from_ref, fill=Category)) + geom_bar(stat='identity') + coord_flip() + xlab("Pathway Name") + ylab(ylabel) + theme(text = element_text(size=15), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10"))
ggplot(ec, aes(x=PathwayName, y=EC_from_ref, fill=Category)) + geom_bar(stat='identity') + coord_flip() + xlab("Pathway Name") + ylab(ylabel) + theme(text = element_text(size=20, face="bold"))
dev.off()

