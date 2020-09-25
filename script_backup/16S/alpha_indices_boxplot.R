# awk -F "\t" '{print $1"\t"$4}' Mapping.txt > group.txt
# Rscript /export/EC1680U/perl/bin/16S/alpha_indices_boxplot.R alpha_indices_table.txt group.txt
# 產生 boxplot.tiff

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

# fn <- 'alpha_indices_table.txt'
# fn2 <- 'group.txt'

fn <- args[1]
fn2 <- args[2]
alph <- read.table(fn, header=TRUE, sep = "\t")
colnames(alph)[1] = 'Sample' # rename column
gp   <- read.table(fn2, sep = "\t") # 不清楚為什麼 header 無法 使用
colnames(gp) <- c("Sample", "Group") # rename column
data <- merge(x = alph, y = gp, by = "Sample", all = TRUE)
tiff(file = "boxplot.tiff", units="in", width=5, height=5, res=300)
par(mfrow=c(1,3))
boxplot(data$chao1 ~ data$Group , col="white", ylab="chao1",  boxwex=0.4 , main="Chao1")
boxplot(data$shannon ~ data$Group , col="white", ylab="shannon",  boxwex=0.4 , main="shannon")
boxplot(data$observed_species ~ data$Group , col="white", ylab="observed_species",  boxwex=0.4 , main="observed_species")
dev.off()
