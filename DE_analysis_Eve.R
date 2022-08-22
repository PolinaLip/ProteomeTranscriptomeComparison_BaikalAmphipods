library(ggplot2)
library(DESeq2)
library(apeglm)
library(ggrepel)
library(dplyr)
library(fgsea)
library(ggfortify) 
library(tibble)

dir_counts <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/Eve/counts' # specify path to your samples
#dir_counts <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/3h/Eve/count'
current_dir <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/Eve'
#current_dir <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/3h/Eve'
files_names <- c('Eve_6C_rep1', 'Eve_6C_rep2', 'Eve_6C_rep3', 'Eve_6C_rep4', 
                 'Eve_24C_rep1', 'Eve_24C_rep2', 'Eve_24C_rep3', 'Eve_24C_rep4') # specify the names of folders with quant.sf data

files_names <- c('Eve_6C_rep1', 'Eve_6C_rep2', 'Eve_6C_rep3', 
                 'Eve_24C_rep1', 'Eve_24C_rep2', 'Eve_24C_rep3', 'Eve_24C_rep4')

countFiles <- list.files(dir_counts, full.names = T)
countFiles
counts <- lapply(countFiles, function(countsFile) {
  f <- file.path(countsFile, 'quant.sf')
  read.table(f, sep="\t", header=1, row.names = 1, 
             stringsAsFactors = F, comment.char = "#")
})
counts <- lapply(counts, function(countsTable) countsTable[, 4, drop=F]) # take NumRead column (the 4th one) -> DESeq performes its own normalization

counts <- do.call(cbind, counts)
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3', '24C_rep4',
                  '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4')
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3', '24C_rep4',
                  '6C_rep1', '6C_rep2', '6C_rep3')
colnames(counts) <- sample_names

coldata <- data.frame(sample=sample_names,
                      temp=factor(c(rep('24C', 4), rep('6C', 4)), 
                                  levels = c('6C', '24C')),
                      row.names=sample_names)
coldata <- data.frame(sample=sample_names,
                      temp=factor(c(rep('24C', 4), rep('6C', 3)), 
                                  levels = c('6C', '24C')),
                      row.names=sample_names)

dds <- DESeqDataSetFromMatrix(countData = round(counts[rowSums(counts) > 10, ]),
                              colData = coldata,
                              design = ~ temp)

dds <- DESeq(dds) # DESeq function requires raw read number -> not TPM. Inside, this function performs its own normalization (median of ratios method) 
resultsNames(dds) # lists the coefficients
deseq_res <- results(dds, name = resultsNames(dds)[2])
t_stat <- deseq_res$stat

vst <- varianceStabilizingTransformation(dds)
plotPCA(vst, intgroup=c("temp")) +
  geom_text_repel(label = vst$sample) +
  theme_light()

res <- lfcShrink(dds, coef="temp_24C_vs_6C", type="apeglm")
head(res)

res <- as.data.frame(res)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=padj < 0.001)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red"))

write.table(rownames_to_column(res, 'contig'),
            file = file.path(current_dir, 'Eve_transcr_24vs6_all_24h_wo__EveB24_4.csv'),
            sep = '\t', quote = F, col.names = T, row.names = F)
# make the table with DE proteins with contig names in the rows, remove .p, take only Ecy?

res[rownames(res) == 'NODE_12622_length_2320_cov_4346.717305_g3383_i2',]

slice_min(res, log2FoldChange, n=5)

res_sign <- subset(res, padj < 0.05)
res_sign <- subset(res_sign, abs(log2FoldChange) > 3)
res_sign_up <- subset(res_sign, log2FoldChange > 0)
res_sign_down <- subset(res_sign, log2FoldChange < 0)

# add annotation 
annot <- read.csv(file.path('~/labeglo2/proteome_transcr_comparision/', 
                            paste0('contigs_whole_annot_Eve','.csv')), sep = '\t') # 24h
counts$annotation <- annot[match(rownames(counts), annot$contig),]$annot
counts_wo_annot <- select(counts, -annotation)
boxplot(counts_wo_annot[rownames(counts_wo_annot) == 'NODE_2787_length_4815_cov_176.045741_g620_i1',])

# pca for samples 

counts_t <- t(counts)
pca_res <- prcomp(counts_t)
set.seed(666)
autoplot(pca_res, data=coldata, colour='temp') +
  geom_text_repel(label = coldata$sample, size = 3) +
  theme_light()
ggsave(file.path(current_dir, 'pca_eve_wo_wrong_sample.png'))

