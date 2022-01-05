library(ggplot2)
library(DESeq2)
library(apeglm)
library(ggrepel)
library(dplyr)
library(fgsea)

dir <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/Gla/counts' # specify path to your samples
dir <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/3h/Gla/count/'
current_dir <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/Gla'
current_dir <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/3h/Gla/'

countFiles <- list.files(dir, full.names = T)
countFiles
counts <- lapply(countFiles, function(countsFile) {
  f <- file.path(countsFile, 'quant.sf')
  read.table(f, sep="\t", header=1, row.names = 1, 
             stringsAsFactors = F, comment.char = "#")
})
counts <- lapply(counts, function(countsTable) countsTable[, 4, drop=F])

counts <- do.call(cbind, counts)
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3', 
                  '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4')
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3', '24C_rep4',
                  '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4')
colnames(counts) <- sample_names

coldata <- data.frame(sample=sample_names,
                      temp=factor(c(rep('24C', 3), rep('6C', 4)), 
                                  levels = c('6C', '24C')),
                      row.names=sample_names)
coldata <- data.frame(sample=sample_names,
                      temp=factor(c(rep('24C', 4), rep('6C', 4)), 
                                  levels = c('6C', '24C')),
                      row.names=sample_names)

dds <- DESeqDataSetFromMatrix(countData = round(counts[rowSums(counts) > 10, ]),
                              colData = coldata,
                              design= ~ temp)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
deseq_res <- results(dds, name = resultsNames(dds)[2])
t_stat <- deseq_res$stat

vst <- varianceStabilizingTransformation(dds)
plotPCA(vst, intgroup=c("temp")) + theme_light()

res <- lfcShrink(dds, coef="temp_24C_vs_6C", type="apeglm")
head(res)

res <- as.data.frame(res)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=padj < 0.05)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red"))

write.table(rownames_to_column(res, 'contig'),
            file = file.path(current_dir, 'Gla_transcr_24vs6_all_3h.csv'),
            sep = '\t', quote = F, col.names = T, row.names = F)
# make the table with DE proteins with contig names in the rows, remove .p, take only Ecy?

res[rownames(res) == 'NODE_12622_length_2320_cov_4346.717305_g3383_i2',]

slice_min(res, log2FoldChange, n=5)


