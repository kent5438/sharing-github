x <- read.delim("FullTable_iso.txt", header=T, sep="\t", check.names=F)
y <- read.delim("Gene_list.txt", header=T, sep="\t", check.names=F)
xy <- merge(x, y, by="Description", all.x=TRUE)
write.table(xy, file="FullTable.iso.txt", sep="\t", quote=F, na="-", row.names = FALSE)
