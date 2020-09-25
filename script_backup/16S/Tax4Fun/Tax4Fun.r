library(Tax4Fun)
library(reshape2)
library(doBy)
library(ggplot2)
#library(Cairo)
#Sys.setlocale(category='LC_ALL', locale='en_US.UTF-8')
folderReferenceData <- '/export/EC1680U/perl/bin/16S/Tax4Fun/SILVA123'
QIIMESingleData <- importQIIMEData('motu_table.taxID.filter.dense.dedup.txt')

##### Gene
Tax4FunOutput <- Tax4Fun(QIIMESingleData, folderReferenceData, fctProfiling = TRUE, refProfile = 'UProC', shortReadMode = TRUE, normCopyNo = TRUE)
tax4fun_gene <- as.data.frame(t(Tax4FunOutput$Tax4FunProfile))
gene <- rownames(tax4fun_gene)
write.table(cbind(gene, tax4fun_gene), 'tax4fun.gene.txt', row.names = FALSE, sep = '\t', quote = FALSE)

##### Pathway
Tax4FunOutput <- Tax4Fun(QIIMESingleData, folderReferenceData, fctProfiling = FALSE, refProfile = 'UProC', shortReadMode = TRUE, normCopyNo = TRUE)
tax4fun_pathway <- as.data.frame(t(Tax4FunOutput$Tax4FunProfile))
pathway <- rownames(tax4fun_pathway)
write.table(cbind(pathway, tax4fun_pathway), 'tax4fun.pathway.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#為 KEGG 第三層級代謝通路結果添加註釋
kegg_anno <- read.delim('/export/EC1680U/perl/bin/16S/Tax4Fun/pathway.anno.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
tax4fun_pathway$Pathway3_level <- t(data.frame(strsplit(rownames(tax4fun_pathway), ';')))[ ,1]
tax4fun_pathway_anno <- merge(tax4fun_pathway, kegg_anno, by = 'Pathway3_level')
# reformat
cols = colnames(tax4fun_pathway_anno)
newcol=cols[c(2:(length(cols)-5), length(cols)-4, length(cols)-3, length(cols)-2, length(cols)-1, 1, length(cols))]
tax4fun_pathway_anno = tax4fun_pathway_anno[newcol]
write.table(tax4fun_pathway_anno, 'tax4fun.pathway.anno.txt', row.names = FALSE, sep = '\t', quote = FALSE)

# 在 KEGG 第二層級統計求和，並添加樣本分組信息
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
names(group)[1] <- 'variable'
tax4fun_pathway_anno_L2 <- tax4fun_pathway_anno[c(group$variable, 'Pathway2_level')]
write.table(tax4fun_pathway_anno_L2, 'tax4fun.pathway.anno.L2.txt', row.names = FALSE, sep = '\t', quote = FALSE)
tax4fun_pathway_anno_L2_melt <- melt(tax4fun_pathway_anno_L2, id = 'Pathway2_level')
write.table(tax4fun_pathway_anno_L2_melt, 'tax4fun.pathway.anno.L2.melt.txt', row.names = FALSE, sep = '\t', quote = FALSE)
tax4fun_pathway_anno_L2_melt_group <- merge(tax4fun_pathway_anno_L2_melt, group, by = 'variable')
write.table(tax4fun_pathway_anno_L2_melt_group, 'tax4fun.pathway.anno.L2.melt.group.txt', row.names = FALSE, sep = '\t', quote = FALSE)

se <- function(x) sd(x) / (length(x))^0.5
pathway_stat_L2 <- summaryBy(value~Description+Pathway2_level, tax4fun_pathway_anno_L2_melt_group, FUN = c(mean, sd, se))
write.table(pathway_stat_L2, 'pathway.stat.L2.txt', row.names = FALSE, sep = '\t', quote = FALSE)

kegg_anno_2 <- kegg_anno[!duplicated(kegg_anno$Pathway2), ][-c(5:6)]
pathway_stat_L2_kegg <- merge(pathway_stat_L2, kegg_anno_2, by = 'Pathway2_level', all.x = TRUE)
# reformat
cols = colnames(pathway_stat_L2_kegg)
newcol=cols[c(2,3,4,5,6,7,1,8)]
pathway_stat_L2_kegg = pathway_stat_L2_kegg[newcol]
write.table(pathway_stat_L2_kegg, 'pathway.stat.L2.kegg.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#使用 ggplot2 作图
pathway_stat_L2_kegg$value.mean <- 100 * pathway_stat_L2_kegg$value.mean
pathway_stat_L2_kegg$value.sd <- 100 * pathway_stat_L2_kegg$value.sd

pathway_stat_L2_kegg$Description = as.character(pathway_stat_L2_kegg$Description)

pathway2_plot = ggplot(pathway_stat_L2_kegg, aes(Pathway2, value.mean, group = Description, fill = Description)) +
geom_col(position = 'dodge', width = 0.8, colour = 'black', size = 0.05) + # dodge 柱狀圖樣式
#geom_errorbar(aes(ymin = 0, ymax = value.mean + value.sd), size = 0.05, width = .35, position = position_dodge(width = .8)) +  #添加误差线（均值±标准差）
coord_flip() + #横、縱坐標軸反轉
theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent',  color = 'black')) +  #去除默认的背景框
# geom_text(aes(label = value.mean, y = value.mean + value.sd + 0.5), size = 4, position = position_dodge(0.8)) +  #添加显著性标记“*”
labs(x = 'KEGG pathway2', y = 'Relative Abundance (%)')  #坐标轴标题

#ggsave('Kegg_pathway2.png', pathway2_plot, dpi=300, width = 20, height = 20, units = "cm")
ggsave('Kegg_pathway2.pdf', pathway2_plot, width = 8, height = 10)
#cairo_pdf('Kegg_pathway2.pdf')
#plot(pathway2_plot)
#dev.off()
