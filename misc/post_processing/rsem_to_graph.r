library(tximport)
library("DESeq2")
dir <- "~/Documents"

#load files into correct format
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE) #samples.txt must be tab delimited
files <- file.path(dir, paste0(samples$run, ".genes.results"))
names(files) <- paste0("sample", 1:2)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)

#samples.txt
#sample  experiment  run
#o.pumilio	skin	epithelium1
#o.pumilio	skin	epithelium2

sampleTable <- data.frame(condition = factor(rep(c("A", "B"))))
rownames(sampleTable) <- colnames(txi.rsem$counts)

library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi.rsem, colData = sampleTable, design = ~ condition)

#differential expression analysis
dds <- DESeq(ddsTxi)
res <- results(dds)
res
