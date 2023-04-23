library("DESeq2")

readcount = "U251_NC_DARS-AS1SI.gene_id.htseqcount.reverse.out"
InputData_counts <- read.delim(readcount, row.names=1)

countTable <- as.matrix(InputData_counts[, dataCols])
head(countTable)
condition <- factor(conditionVec)
coldata <- data.frame(row.names=colnames(countTable), condition)
coldata

dds <- DESeqDataSetFromMatrix(countData=countTable, colData=coldata, design= ~condition)
dds <- DESeq(dds)
dds$sizeFactor

res <- results(dds)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Additional output annotation info and raw count
resdata_rawdata <- merge(as.data.frame(InputData_counts[, c("Symbol", "Type", dataCols)]), as.data.frame(resdata), by.x="row.names", by.y="Gene", sort=FALSE)
names(resdata_rawdata)[1] <- "Gid"
head(resdata_rawdata)
write.csv(resdata_rawdata, file=paste(c(OutputPrefix, DE_groups, "diffexpr-results.csv"), collapse="_"), row.names = FALSE)
