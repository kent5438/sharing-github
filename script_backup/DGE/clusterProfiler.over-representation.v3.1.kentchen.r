suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(pathview))

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

fn <- args[1]
genome<-args[2]
dir<-args[3]

effect  <- read.table (fn, header = TRUE)

# add species control structure 
if ((genome == "hg19") || (genome == "human")) {
  genome = "human"
  OrgDb = "org.Hs.eg.db" 
  species = "hsa"
}
if ((genome == "mm10") || (genome == "mm9")) {
  genome = "mouse"
  OrgDb = "org.Mm.eg.db"
  species = "mmu"
}
if (genome == "rat") {
  genome = "rat"
  OrgDb = "org.Rn.eg.db"
  species = "rno"
}

gen.entrzid <- bitr(row.names(effect), fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb, drop=FALSE)

### Biological Process
output <- paste0(dir, "/gen.bp.enrichGO.txt");
gen.bp.ego <- enrichGO(gene = gen.entrzid$ENTREZID, OrgDb=OrgDb, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.2, readable = TRUE)
#gen.bp.ego <- enrichGO(gene = gen.entrzid$ENTREZID, OrgDb=OrgDb, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.3, readable = TRUE)
write.table(summary(gen.bp.ego) , file = output, append = FALSE, quote = FALSE, sep = "\t",   eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double", fileEncoding = "")

### dotplot
output = paste0(dir, "/gen.bp.ego.dotplot.png");
png(filename = output, width = 1600, height = 1200)
try(dotplot(gen.bp.ego, font.size = 12), TRUE)
dev.off()

###  plotGOgraph
output = paste0(dir, "/gen.bp.ego.plotGOgraph.pdf");
pdf (file = output)
try(plotGOgraph(gen.bp.ego), TRUE)
dev.off()

### Molecular Function
output = paste0(dir, "/gen.mf.enrichGO.txt");
gen.mf.ego <- enrichGO(gene = gen.entrzid$ENTREZID, OrgDb=OrgDb, ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.2, readable = TRUE)
#gen.mf.ego <- enrichGO(gene = gen.entrzid$ENTREZID, OrgDb=OrgDb, ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.3, readable = TRUE)
write.table(summary(gen.mf.ego) , file = output, append = FALSE, quote = FALSE, sep = "\t",   eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double", fileEncoding = "")

### dotplot
output = paste0(dir, "/gen.mf.ego.dotplot.png");
png(filename = output, width = 1600, height = 1200)
try(dotplot(gen.mf.ego, font.size = 12), TRUE)
dev.off()

###  plotGOgraph
output = paste0(dir, "/gen.mf.ego.plotGOgraph.pdf");
pdf (file = output)
try(plotGOgraph(gen.mf.ego), TRUE)
dev.off()


### Cellular Component
output = paste0(dir, "/gen.cc.enrichGO.txt"); 
gen.cc.ego <- enrichGO(gene = gen.entrzid$ENTREZID, OrgDb=OrgDb, ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.2, readable = TRUE)
#gen.cc.ego <- enrichGO(gene = gen.entrzid$ENTREZID, OrgDb=OrgDb, ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.3, readable = TRUE)
write.table(summary(gen.cc.ego) , file = output, append = FALSE, quote = FALSE, sep = "\t",   eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double", fileEncoding = "")

### dotplot
output = paste0(dir, "/gen.cc.ego.dotplot.png"); 
png(filename = output, width = 1600, height = 1200)
try(dotplot(gen.cc.ego, font.size = 12), TRUE)
dev.off()

###  plotGOgraph
output = paste0(dir, "/gen.cc.ego.plotGOgraph.pdf"); 
pdf (file = output)
try(plotGOgraph(gen.cc.ego), TRUE)
dev.off()

### KEGG section

output = paste0(dir, "/gen.enrichKEGG.txt")
gen.kegg <- enrichKEGG(gen.entrzid$ENTREZID, organism = species, pvalueCutoff  = 0.05, qvalueCutoff  = 0.2, pAdjustMethod="BH")
#gen.kegg <- enrichKEGG(gen.entrzid$ENTREZID, organism = species, pvalueCutoff  = 0.3, pAdjustMethod="BH")
write.table(summary(gen.kegg) , file = output, append = FALSE, quote = FALSE, sep = "\t",   eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = "double", fileEncoding = "")

### KEGG PathView
target.workdir<-paste0(dir, "/keggmap")
dir.create(target.workdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
setwd(target.workdir)

log2FC<-as.vector(log2(effect$PostFC))
# only keep the unique values in SYMBOL column
names(log2FC)<-gen.entrzid[!duplicated(gen.entrzid[,c('SYMBOL')]),]$ENTREZID

gen.kegg.pathway<-as.vector(summary(gen.kegg)$ID)
if (length(rownames(summary(gen.kegg)))!=0){
  for (x in 1:length(gen.kegg.pathway)){
    #temp <- pathview(gene.data  = log2FC, pathway.id = gen.kegg.pathway[x], species = species, limit = list(gene=max(abs(log2FC)), cpd=1))
	try(pathview(gene.data  = log2FC, pathway.id = gen.kegg.pathway[x], species = species, limit = list(gene=max(abs(log2FC)), cpd=1)), TRUE)
  }
}
