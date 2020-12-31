library(edgeR)
library(limma)
library(Glimma)
DE_results <- read.delim("Trans.raw.counts.matrix.sample1_vs_sample2.edgeR.DE_results", row.names=1, stringsAsFactors=FALSE, check.names=F)
DE_results <- DE_results[order(row.names(DE_results)),]
rawcounts <- read.table("Trans.raw.counts.matrix", header=T, row.names=1, com='', check.names=F)
rawcounts <- rawcounts[order(row.names(rawcounts)),]
samples <- read.delim("sample.txt", header=FALSE, stringsAsFactors=FALSE, com='#', check.names=F)
colnames(samples) <- c("Group","Sample")
#rawcounts <- rawcounts[,samples$Sample]
rawcounts[rawcounts ==0] <- 0.000000001
rnaseqMatrix = rawcounts
#rnaseqMatrix = round(rawcounts)
#rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=2,]
conditions = factor(samples$Group)
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study$common.dispersion=0.1
exp_study = estimateCommonDisp(exp_study)
exp_study$common.dispersion=0.1
exp_study = estimateTagwiseDisp(exp_study)
cpm_table <- cpm(exp_study)

cpm_table <- as.data.frame(cpm_table)
cpm_table = cpm_table[rownames(DE_results),]


### Optional, Merge the DE_results and cpm_table
cpm.DE <- cbind(cpm_table, DE_results)
write.csv(cpm.DE, file="cpm.DE_results.csv")

### Set significance level for different colors on plot
DE_results$GeneID <- rownames(DE_results)
DE_results$logCPM <- DE_results$logCPM
DE_results$Sig <- as.numeric(DE_results$FDR < 0.05)


### Make the Glimma MA-plot
glMDPlot(DE_results, counts=cpm_table, samples=samples$Sample,
		 anno=cpm.DE, groups=samples$Group,
         xval="logCPM", yval="logFC",
         display.columns=colnames(DE_results),
         folder="glimma_MA", html="MA-plot",
         search.by="GeneID", status=DE_results$Sig, cols=c("black","red"), launch=FALSE)

### Make the Glimma Volcano plot
glXYPlot(x=DE_results$logFC, y=-log10(DE_results$FDR), counts=cpm_table, samples=samples$Sample,
         anno=DE_results, groups=samples$Group,
         xlab="log2(FC)", ylab="-log10(FDR)",
         display.columns=colnames(DE_results),
         folder="glimma_Volcano", html="Volcano-plot",
         search.by="GeneID", status=DE_results$Sig, cols=c("black","red"), launch=FALSE)
