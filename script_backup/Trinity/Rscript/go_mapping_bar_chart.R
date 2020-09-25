# Rscript /export/EC1680U/Rscript/Pacbio/kegg_horizontal_bar_chart.R 5_KEGG/PB17012_NAI5_ec2kegg.txt 5_KEGG ead
# Author: kentchen
# Date: 2017-06-12

library(ggplot2)

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
df <- args[1]
dir <- args[2]

ec <- read.table(df, header=T, sep="\t", check.names=F)
ec <- ec[(ec$Count) > 0,]
ec$Function <- factor(ec$Function, as.character(ec$Function))
output = paste0(dir, "/GO_barchart.png")

png(output, width = 25, height = 20, units = "in", res = 300)
ggplot(ec, aes(x=Function, y=Count, fill=Class)) + geom_bar(stat="identity") + coord_flip() + theme(text = element_text(size=20, face="bold"))
dev.off()
