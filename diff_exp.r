
diff_exp <- function(filename, coldata_filename) {
  #Differential expression analysis
  
  library(EDASeq)
  library(DESeq2)
  
  #Reading the table from file
  filename <- "out"
  coldata_filename <- "coldata.tsv"
  mirna <- read.table(filename, header=TRUE, sep="\t", row.names=1)
  coldata <- read.table(coldata_filename, header=TRUE, sep="\t", row.names=1)
  
  #Removing all the useless column
  mirna <- mirna[,-c(1, 2, 3, 4, 5)]
  
  #Filtering out low expressed miRNA using classical approach
  filter <- apply(mirna, 1, function(x) length(x[which(x>10)])>0)
  filtered_mirna <- as.matrix(mirna)[filter,]
  
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
  return(res);
}