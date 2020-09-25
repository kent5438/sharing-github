# Rscript /export/EC1680U/Rscript/Pacbio/cog_horizontal_bar_chart.R 4_COG/results/func_stats.txt 4_COG/results
# Author: kentchen
# Date: 2017-06-01

library(ggplot2)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
fn<-args[1]
dir<-args[2]

df <- read.table(fn, header=F, sep="\t", check.names=F)
df$V3 = df$V3 / sum(df$V3) * 100
colnames(df) <- c("Class", "Function", "Count")
df <- df[(df$Count) > 0,]
#df$Count <- log2(df$Count)
df <- df[order(df$Count),]
df$Function <- factor(df$Function, as.character(df$Function))

output = paste0(dir, "/cog_barchart.png");
png(output, units="in" , width=25, height=20, res=150)
ggplot(df, aes(x=Function, y=Count, fill=Class)) + geom_bar(stat = 'identity') + coord_flip() + theme(text = element_text(size=20, face="bold")) + xlab("Function") + ylab("Count (%)")
dev.off()
