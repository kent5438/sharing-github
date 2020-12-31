#! /usr/bin/env R

library(ggplot2)
df <- read.table("PCoA_weighted_nmds.2d.txt", header=T)
a <- ggplot(df, aes(x=NMDS1, y=NMDS2, col=Description))+geom_point() +stat_ellipse() +theme_bw() +labs(title = "NMDS Plot")
ggsave("PCoA_weighted_nmds.2d.noSampleName.png", plot=a)
