library(ggplot2)
library(ggfortify)


args <- commandArgs(trailingOnly = TRUE)
argv1 <- args[1]
pca_input <- paste0("meta16S_OUTPUT/otu_table.forPCA.txt")
pca_output <- paste0("meta16S_OUTPUT/PCA.pdf")


df <- read.table(file = pca_input, header = T, sep = "\t", check.names = F)
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
ggsave( filename = pca_output, plot = pca_plot, width = 7, height = 7, units = "in")

