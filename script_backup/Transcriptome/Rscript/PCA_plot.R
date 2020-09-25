library(ggplot2)
library(ggfortify)


args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
out <- args[2]

df <- read.table(file = fn, header = T, sep = "\t", check.names = F)
df <- df[,c(2:ncol(df))]
df <- as.data.frame(t(df))
df <- scale(df)
df[is.na(df)] <- 0
df <- df[,apply(df, 2, var, na.rm=TRUE) != 0]
df_pca <- prcomp(df)
df_out <- as.data.frame(df_pca$x)
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
df_pca_out <- prcomp(df_out)
pca_plot <- autoplot(df_pca_out, label=TRUE, label.vjust=2, label.size=2, type="obs", ggplot2.xlab=percentage[1], ggplot2.ylab=percentage[2])
#pca_plot <- autoplot(df_pca_out, label=TRUE, label.vjust=2, label.size=2, type="obs")
ggsave( filename = out, plot = pca_plot, width = 7, height = 7, units = "in")

