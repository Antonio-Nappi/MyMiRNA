#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(EDASeq)
library(DESeq2)
library(gplots)

#Reading the table from file
filename <- "mirna_out"
coldata_filename <- "coldata.tsv"
#filename <- args[1]
#coldata_filename <- args[2]

#pvalue_filter <- args[3]
pvalue_filter = 0.98
if(pvalue_filter == "none") pvalue_filter <- -1
#padj_filter <- args[4]
padj_filter = "none"
if(padj_filter == "none") padj_filter <- -1
#logfold_filter <- args[5]
logfold_filter = "none"
if(logfold_filter == "none") logfold_filter <- -1000000
#index <- args[6]
index = "Prova"
mirna <- read.table(filename, header=TRUE, sep="\t", row.names=1)
coldata <- read.table(coldata_filename, header=TRUE, sep="\t", row.names=1)

#Removing all the useless column
mirna <- mirna[,-c(1, 2, 3, 4, 5)]

#Filtering out low expressed miRNA using classical approach
filter <- apply(mirna, 1, function(x) length(x[which(x>10)])>0)
filtered_mirna <- as.matrix(mirna)[filter,]
names = as.matrix(row.names(filtered_mirna))
#write(names, file=paste(index,"/mirna/mirna.names", sep=""), ncolumns = 1)

#Normalization
normalized_mirna <- betweenLaneNormalization(filtered_mirna, which="upper")

#Differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = normalized_mirna, 
                              colData = coldata,
                              design = ~Condition)

dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)
p = res[res$pvalue >= pvalue_filter & res$padj >= padj_filter & res$log2FoldChange > logfold_filter, ]

#jpeg("./file.jpg", width = 1200, height = 1200)
#heatmap.2(as.matrix(p),scale = "column", col = heat.colors(256),main = "Test",dendrogram = "column",margins = c(8,10),cexCol = 1)
#dev.off()
test <-t(normalized_mirna)
#barplot(as.matrix(test)/sum(test)*100,legend.text=c("% of miRNA in the wild-type sample","% of miRNA in the tumoral sample"),main="Barplot",col=c('red','green'),xlab="miRNA",ylab="number of reads",beside=FALSE)
