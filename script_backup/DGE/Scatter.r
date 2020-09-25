# Rscript Scatter.r output/ebseq/AM1-D7,AM3-D7-AM1-2D-1.genMat.result.1 ./
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
fn<-args[1]
dir<-args[2]

data <- read.table(fn, header = TRUE, sep = "\t")
# colnames(data)
C1Mean <-data$C1Mean
C2Mean <-data$C2Mean
PValue <-data$PPDE

col = rep('black', length(C1Mean))
col[PValue > 0.95] <- 'blue'


output = paste0(dir, "/Scatter.png");
png(filename = output, width = 1600, height = 1200)
par(col="black")
plot (log2(C1Mean+1), log2(C2Mean+1), col=col, main="Scatter plot of C1Mean and C2Mean", xlab="log2(C1Mean)", ylab="log2(C2Mean)")
# add two red straight lines
abline(1,1, col = 2)
abline(-1,1, col = 2)
# add regression line
abline(lm(log2(C1Mean+1)~log2(C2Mean+1)), col=3)
dev.off()

