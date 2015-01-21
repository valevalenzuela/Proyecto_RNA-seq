> source("http://bioconductor.org/biocLite.R") 
> biocLite("DESeq2") 

> library(DESeq2)

> data <- read.table("Desktop/Proyecto_RNA-seq/data/Cgigas-HS-count.txt", header = T, sep = "\t")
> rownames(data) <- data$Feature
> data <- data[,-1]

comando para comparar dos set de datos (volcano)
> deseq2.colData <- data.frame(condition=factor(c(rep("Treated", 3), rep("Control", 3))), type=factor(rep("heat shock", 6)))
> rownames(deseq2.colData) <- colnames(data)
> deseq2.dds <- DESeqDataSetFromMatrix(countData = data, colData = deseq2.colData, design = ~ condition)
> deseq2.dds <- DESeq(deseq2.dds)

> deseq2.res <- results(deseq2.dds)
> deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]
> head(deseq2.res)

> dim(deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ])
> plot(tmp$baseMean, tmp$log2FoldChange, pch=20, cex=0.45, ylim=c(-3, 3), log="x", col="blue", main="DEG Virus Exposure  (pval <= 0.05)", xlab="mean of normalized counts", ylab="Log2 Fold Change")
> tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ]
> points(tmp.sig$baseMean, tmp.sig$log2FoldChange, pch=20, cex=0.45, col="red")
> abline(h=c(-1,1), col="black")
> write.table(tmp.sig, "Desktop/Proyecto_RNA-seq/output/Cgigas.tab", row.names = T)