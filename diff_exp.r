#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(EDASeq)
library(DESeq2)
library(gplots)
library(RColorBrewer)

#Reading the table from file
#filename <- "mirna_out"
#coldata_filename <- "coldata.tsv"
filename <- args[1]
coldata_filename <- args[2]

pvalue_filter <- args[3]
#pvalue_filter = "0.05"
if(pvalue_filter == "none") pvalue_filter <- 2
padj_filter <- args[4]
#padj_filter = "none"
if(padj_filter == "none") padj_filter <- 2
logfold_filter <- args[5]
#logfold_filter = "none"
if(logfold_filter == "none") logfold_filter <- -1000000
index <- args[6]
index = "prova1"

mirna <- read.table(filename, header=TRUE, sep="\t", row.names=1)
coldata <- read.table(coldata_filename, header=TRUE, sep="\t", row.names=1)
palette <- brewer.pal(8, 'Set2')

#Removing all the useless column
mirna <- mirna[,-c(1, 2, 3, 4, 5)]

#Filtering out low expressed miRNA using classical approach
filter <- apply(mirna, 1, function(x) length(x[which(x>10)])>0)
filtered_mirna <- as.matrix(mirna)[filter,]
names = as.matrix(row.names(filtered_mirna))
write(names, file=paste(index,"/mirna/mirna.names", sep=""), ncolumns = 1)

#Normalization
normalized_mirna <- betweenLaneNormalization(filtered_mirna, which="upper")

#BOX PLOT
jpeg(paste("gui/src/assets/",index,"/boxplot.jpeg", sep=""), width = 1280, height = 720)
plotRLE(normalized_mirna,
        outline=FALSE, las=3, col=palette,
        ylab="Relative Log Expression",
        cex.axis=0.7,
        cex.lab=1)
dev.off()

#Differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = normalized_mirna, 
                              colData = coldata,
                              design = ~Condition)

dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)
filtered_res = res[res$pvalue <= pvalue_filter & res$padj <= padj_filter & res$log2FoldChange > logfold_filter, ]
print(dim(res))
print(dim(filtered_res))
#SPREAD PLOTS
jpeg(paste("gui/src/assets/",index,"/spread_all.jpeg", sep=""), width = 1280, height = 720)
DESeq2::plotMA(res, ylim=c(-5, 5))
dev.off()

jpeg(paste("gui/src/assets/",index,"/spread_filtered.jpeg", sep=""), width = 1280, height = 720)
DESeq2::plotMA(filtered_res, ylim=c(-5, 5))
dev.off()

# HISTOGRAMS
jpeg(paste("gui/src/assets/",index,"/hist_p_all.jpeg", sep=""), width = 1280, height = 720)
hist(res$pvalue)
dev.off()

jpeg(paste("gui/src/assets/",index,"/hist_p_filtered.jpeg", sep=""), width = 1280, height = 720)
hist(filtered_res$pvalue)
dev.off()

jpeg(paste("gui/src/assets/",index,"/hist_padj_all.jpeg", sep=""), width = 1280, height = 720)
hist(res$padj)
dev.off()

jpeg(paste("gui/src/assets/",index,"/hist_padj_filtered.jpeg", sep=""), width = 1280, height = 720)
hist(filtered_res$padj)
dev.off()

# HEATMAP
jpeg(paste("gui/src/assets/",index,"/heatmap.jpeg", sep=""), width = 1280, height = 720)
heatmap.2(as.matrix(filtered_res),scale = "column", col = heat.colors(256),main = "Test",dendrogram = "column",margins = c(8,10),cexCol = 1)
dev.off()

# BAR PLOT
jpeg(paste("gui/src/assets/",index,"/barplot.jpeg", sep=""), width = 1280, height = 720)
test <-t(normalized_mirna[row.names(filtered_res), ])
barplot(as.matrix(test)/sum(test)*100,
        legend.text=c("% of miRNA in the wild-type sample",
                      "% of miRNA in the tumoral sample"),main="Barplot",
        col=palette ,xlab="miRNA",ylab=" % number of reads",beside=FALSE,
        cex.names=0.5)
dev.off()

# VOLCANO PLOT
jpeg(paste("gui/src/assets/",index,"/volcano.jpeg", sep=""), width = 1280, height = 720)
plot(res$log2FoldChange, -log10(res$pvalue), 
     xlab = "log2FC", ylab = "-log10(p)", pch=20, col="grey")
if (padj_filter == 2) padj_filter <- 0.05
de <- which(res$padj <= padj_filter)
points(res[de, "log2FoldChange"],
       -log10(res[de, "pvalue"]),
       pch=20, col=3)
legend("topright", legend=c("All miRNAs",
                       paste("miRNAs with an adjusted p-value under", padj_filter)),
       col=c("grey", "green"), cex=0.8, pch=c(1,1))
dev.off()