# DE analysis for the each group with missing values (without NAs, with 3 NAs 
# (3 TMT batch), 8 NAs (1 TMT batch), 8 NAs (2 TMT batch), 11 NAs (1 and 3 TMT 
# batches), 11 NAs (1 and 2 TMT batches))

library(tidyverse)
library(proteus)
library(limma) 
library(edgeR) 
library(sva) # from Bioconductor
library(psych) # library(devtools) -> install_version('mnormt', '1.5-7') -> install.packages("psych", repos = "https://personality-project.org/r/", type="source")
library(qvalue)
library(pheatmap)
library(viridis)
library(ggplot2)
library(ggfortify)
library(ggpol)
library(imp4p)

### 1. Upload metafile
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  meta <- meta[grepl(species_name, meta$condition),]
  return(meta)
}
path2meta <-
  'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/Metadata_Proteus.tsv'

meta <- meta_upload(path2meta, 'Eve')

### 2. Data uploading (proteinGroups file from MaxQuant)
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_eve/'
proteinGroups_file <- 'proteinGroups_wo_cont_more2pept.txt' # take file with proteinGroups with 2 or more peptides quantified
dat_init <- read.csv(file.path(dir, proteinGroups_file), sep = '\t', header = T, 
                     check.names = F) 
pep_annot <- read.csv(file.path(dir, 'annot_protein_groups_eve.csv'), sep = '\t',
                      header = T)

select_data <- function(meta_data, proteinGroups_data){
  dat <- proteinGroups_data[,c('Protein IDs', meta_data$measure)]
  rownames(dat) <- dat$`Protein IDs`
  dat <- dat[-1]
  colnames(dat) <- meta$sample
  dat[dat == 0] <- NA
  return(dat)
}
dat <- select_data(meta, dat_init)

### 3. Get rid of the outliers and samples, which you do not want to analize
dat <- dat[,!grepl('VrK1_390_3|VrK1_390_4', colnames(dat))]
dat <- dat[,!grepl('Eve_pool', colnames(dat))]
meta <- subset(meta, sample != 'VrK1_390_3' & sample != 'VrK1_390_4' &
                 !grepl('Eve_pool', sample)
)

### 4. Remove raws with all NAs in at least one of the condition
cond <- as.factor(meta$condition)
ixs <- order(cond)
dat <- dat[, ixs]
cond <- cond[ixs]

for (level in levels(cond)) {
  col_ixs <- cond == level
  dat <- dat[rowSums(is.na(dat[,col_ixs])) < sum(col_ixs),]
}

data_raw <- dat

### 5. Visualize the raw data
boxplot(log2(data_raw), col = c(rep('blue', 8), rep('red', 5), rep('green', 6)), 
        notch = TRUE, main = 'RAW data: 24C (blue), 6C_after (red), 6C_before (green)',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_raw), col = c(rep('blue', 8), rep('red', 5), rep('green', 6)), 
              main = 'Raw data', legend=F) # Eve, resolution -> 800/553

### 6. Sample loading normalization

exp1_raw <- data_raw[,subset(meta, experiment == 6)$sample]
exp2_raw <- data_raw[,subset(meta, experiment == 7)$sample]
exp3_raw <- data_raw[,subset(meta, experiment == 8)$sample]
target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T), 
               colSums(exp3_raw, na.rm = T)),
               na.rm = T)
norm_facs <- target / colSums(exp1_raw, na.rm = T)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw, na.rm = T)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3_raw, na.rm = T)
exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")

data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl)

boxplot(log2(data_sl), col = c(rep(c('red', 'blue'), each = 8), rep('green', 3)), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nExp1 (red), Exp2 (blue), Exp3 (green)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_sl), col = c(rep(c('red', 'blue'), each = 8), rep('green', 3)),
              main = "SL normalization", legend = F) # Eve, resolution -> 800/553

### 7. Vizualize normalized data
## 7.1 MDS-plot with groups colored by TMT batch
plotMDS(log2(data_sl), col = c(rep(c('red', 'blue'), each = 8), rep('green', 3)), 
        main = "SL clusters grouped by TMT experiment") # Eve

## 7.2 MDS-plot with groups colored by condition:
png(file.path(dir, 'Eve_MDS_wIMBR_DreamAIimput.png'),
    width = 1200, height = 1000, pointsize = 28)
col_vec <- ifelse(grepl('pool', colnames(data_sl)), 'green',
                  ifelse(grepl('VrK1', colnames(data_sl)), 'orange',
                         ifelse(grepl('Vr1', colnames(data_sl)), 'blue', 'red')))

plotMDS(log2(data_sl), col = col_vec)
dev.off()

## 7.3 PCA-plot
data_sl_t <- na.omit(data_sl) %>% t %>% as.data.frame
pca_res <- prcomp(data_sl_t)
meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
autoplot(pca_res, data=meta2, colour='condition') + theme_light()
ggsave(filename = file.path(dir, 'pca_proteinGroups_woImput_woPools.png'))

### 8. Multiply the normalized data by 10^7 to make it possible to perform DE analysis
data_sl_mult <- data_sl * 10000000 # for PSM-normalized data
data_irs <- data_sl_mult # for PSM-normalized data

#########
# EdgeR 
#########
### 9. Separate the data on the subset with the same NA pattern and perform EdgeR analysis

# here are two functions to perform de analysis
edger_analysis <- function(datatode, group, time_group, cond_vector, 
                           protein_annotation){
  
  edge_object <- DGEList(counts = datatode, group = group)
  design <- model.matrix(~0 + group + time_group)# + group_outlier)
  edge_object <- estimateDisp(edge_object, design)
  plotBCV(edge_object, main = "Biological variation SL only (no IRS)")
  edge_res <- exactTest(edge_object, pair = cond_vector)
  
  print(summary(decideTestsDGE(edge_res)))
  
  tt_edge_res <- topTags(edge_res, n = Inf, sort.by = "none")
  tt_edge_res <- tt_edge_res$table
  tt_edge_res_sign <- subset(tt_edge_res, PValue < 0.05)
  
  proteins <- rownames(tt_edge_res)
  samples2annot <- match(proteins, protein_annotation$protein_group)
  sum(is.na(samples2annot))
  
  SignProteinInfo <- data.frame(protein = proteins,
                                geneSymbol = protein_annotation$annotation[samples2annot],
                                logFC = tt_edge_res$logFC,
                                pvalue = tt_edge_res$PValue,
                                FDR = tt_edge_res$FDR)
  
  SignProteinInfo <- SignProteinInfo[order(SignProteinInfo$FDR),]
  return(SignProteinInfo)
}

edger_analysis_wrap <- function(norm_data, protein_annotation){
  data2DE <- norm_data[,grepl('VrK2|Vr1|VrK1', colnames(norm_data))] # 24C vs 6C
  group <- ifelse(grepl('VrK1|VrK2', colnames(norm_data)), '6C', '24C')
  group_time <- as.factor(ifelse(grepl('VrK1|CK1|LK1', colnames(data2DE)), 25, 26))
  return(edger_analysis(data2DE, group, group_time, c('6C', '24C'), 
                        protein_annotation))
}

# wo NAs
data_wo_na <- na.omit(data_irs)
SignProteinInfo_0 <- edger_analysis_wrap(data_wo_na, pep_annot)

# 3 NAs/row
data_3 <- data_irs[rowSums(is.na(data_irs)) == 3,]
data_3 <- data_3[,colSums(is.na(data_3))<nrow(data_3)]
SignProteinInfo_3 <- edger_analysis_wrap(data_3, pep_annot)

# 8 NAs/row (the 1st + 3rd TMT batches)
data_8_1 <- data_irs[rowSums(is.na(data_irs)) == 8,]
data_8_1 <- data_8_1[,grepl('_1$|_2$|_3$|Vr1_390_7|Vr1_390_8|VrK1_390_8', 
                          colnames(data_8_1))]
data_8_1 <- na.omit(data_8_1)
SignProteinInfo_8_1 <- edger_analysis_wrap(data_8_1, pep_annot)

# 8 NAs/row (the 2nd + 3rd TMT batches)
data_8_2 <- data_irs[rowSums(is.na(data_irs)) == 8,]
data_8_2 <- data_8_2[,grepl('_4$|_5$|_6$|VrK1_390_7|Vr1_390_7|Vr1_390_8|VrK1_390_8', 
                            colnames(data_8_2))]
data_8_2 <- na.omit(data_8_2)
SignProteinInfo_8_2 <- edger_analysis_wrap(data_8_2, pep_annot)

# 11 NAs/row (the 1st TMT batch)
data_11_1 <- data_irs[rowSums(is.na(data_irs)) == 11,]
data_11_1 <- data_11_1[,!grepl('Vr1_390_7|Vr1_390_8|VrK1_390_8', 
                             colnames(data_11_1))]
data_11_1 <- data_11_1[,grepl('_1$|_2$|_3$', 
                            colnames(data_11_1))]
data_11_1 <- na.omit(data_11_1)
SignProteinInfo_11_1 <- edger_analysis_wrap(data_11_1, pep_annot)

# 11 NAs/row (the 2nd TMT batch)
data_11_2 <- data_irs[rowSums(is.na(data_irs)) == 11,]
data_11_2 <- data_11_2[,!grepl('Vr1_390_7|Vr1_390_8|VrK1_390_8', 
                             colnames(data_11_2))]
data_11_2 <- data_11_2[,grepl('_4$|_5$|_6$|VrK1_390_7', 
                            colnames(data_11_2))]
data_11_2 <- na.omit(data_11_2)
SignProteinInfo_11_2 <- edger_analysis_wrap(data_11_2, pep_annot)

### 10. Combine all subsets in one and calculate the adjusted p-values
SPI_all <- rbind(SignProteinInfo_0, SignProteinInfo_3, SignProteinInfo_8_1,
                 SignProteinInfo_8_2, SignProteinInfo_11_1, SignProteinInfo_11_2)

SPI_all$FDR_recalc <- p.adjust(SPI_all$pvalue, method = "fdr")

SPI_all_sign <- subset(SPI_all, FDR_recalc < 0.05)

nrow(subset(SPI_all_sign, logFC > 0)) # UP
nrow(subset(SPI_all_sign, logFC < 0)) # DOWN

########################################
# 11. Draw boxplots for the DE proteins
########################################
res_sign_fdr_up <- subset(SPI_all_sign, logFC > 0)
res_sign_fdr_down <- subset(SPI_all_sign, logFC < 0)
dat_mult_stand_fdr <- apply(data_sl_mult, 1, scale) # if take multiplied and scaled data to vizualize
dat_mult_stand_fdr <- t(dat_mult_stand_fdr) # --||--
#dat_mult_stand_fdr <- data_sl # if take "naked" data after normalization
colnames(dat_mult_stand_fdr) <- colnames(data_sl_mult)
dat_sign_fdr <- dat_mult_stand_fdr[row.names(dat_mult_stand_fdr) %in% 
                                     SPI_all_sign$protein,]

dat_sign_fdr_up <- dat_mult_stand_fdr[row.names(dat_mult_stand_fdr) %in% 
                                        res_sign_fdr_up$protein,]
dat_sign_fdr_down <- dat_mult_stand_fdr[row.names(dat_mult_stand_fdr) %in% 
                                          res_sign_fdr_down$protein,]
prot_annot_up <- match(rownames(dat_sign_fdr_up), pep_annot$protein_group)
prot_annot_down <- match(rownames(dat_sign_fdr_down), pep_annot$protein_group)

gene_names_up <- pep_annot$annotation[prot_annot_up]
gene_names_down <- pep_annot$annotation[prot_annot_down]
gene_names_up <- as.vector(gene_names_up)
gene_names_down <- as.vector(gene_names_down)

dat_sign_fdr_up <- as.data.frame(dat_sign_fdr_up)
dat_sign_fdr_down <- as.data.frame(dat_sign_fdr_down)

dat_sign_fdr_up$protein_name <- gene_names_up
dat_sign_fdr_up$protein_name <- sub('PREDICTED: |LOW QUALITY PROTEIN: ', '', dat_sign_fdr_up$protein_name)
dat_sign_fdr_down$protein_name <- gene_names_down

dat_sign_fdr_up$contig <- rownames(dat_sign_fdr_up)
dat_sign_fdr_up_long <- pivot_longer(data = dat_sign_fdr_up, cols = !protein_name &
                                       !contig)
dat_sign_fdr_up_long$Condition <- ifelse(grepl('Vr1', dat_sign_fdr_up_long$name), 
                                         '24C', '6C')

var_width <- 25
dat_sign_fdr_up_long <- mutate(dat_sign_fdr_up_long, 
                               pretty_varname = str_wrap(dat_sign_fdr_up_long$protein_name, 
                                                         width = var_width))
contig_to_protein <- aggregate(pretty_varname ~ contig,
                               dat_sign_fdr_up_long, head, n = 1)
contig_to_protein2 <- contig_to_protein$pretty_varname
names(contig_to_protein2) <- contig_to_protein$contig

ggplot(dat_sign_fdr_up_long, aes(x = Condition, y = value)) +
  geom_boxplot(fill = 'white', outlier.alpha = 0) +
  geom_jitter(aes(color = Condition), size = 0.8) +
  facet_wrap(~contig,
             labeller = labeller(contig = function(x) contig_to_protein2[x]), 
             ncol = 9) +
  ylab('Scaled intensities') +
  xlab('') +
  theme_bw() +
  theme(strip.text = element_text(size=6),
        legend.position = "none")

# STOP, think
ggsave(filename = file.path(dir, 'DEup_combinedEdgeR_Scaled.png'),
       width = 12, height = 6)
ggsave(filename = file.path(dir, 'DEup_combinedEdgeR_woScaling.png'),
       width = 12, height = 6)

# and now for down-regulated:

dat_sign_fdr_down$protein_name <- gene_names_down
dat_sign_fdr_down$protein_name <- sub('PREDICTED: |LOW QUALITY PROTEIN: ', '', 
                                      dat_sign_fdr_down$protein_name)

dat_sign_fdr_down$contig <- rownames(dat_sign_fdr_down)
dat_sign_fdr_down_long <- pivot_longer(data = dat_sign_fdr_down, 
                                       cols = !protein_name &
                                       !contig)
dat_sign_fdr_down_long$Condition <- ifelse(grepl('Vr1', dat_sign_fdr_down_long$name), 
                                         '24C', '6C')

var_width <- 25
dat_sign_fdr_down_long <- mutate(dat_sign_fdr_down_long, 
                          pretty_varname = str_wrap(dat_sign_fdr_down_long$protein_name, 
                                                         width = var_width))
contig_to_protein_down <- aggregate(pretty_varname ~ contig,
                               dat_sign_fdr_down_long, head, n = 1)
contig_to_protein_down2 <- contig_to_protein_down$pretty_varname
names(contig_to_protein_down2) <- contig_to_protein_down$contig

ggplot(dat_sign_fdr_down_long, aes(x = Condition, y = value)) +
  geom_boxplot(fill = 'white', outlier.alpha = 0) +
  geom_jitter(aes(color = Condition), size = 0.8) +
  facet_wrap(~contig,
             labeller = labeller(contig = function(x) contig_to_protein_down2[x]), 
             ncol = 7) +
  ylab('Intensities') +
  xlab('') +
  theme_bw() +
  theme(strip.text = element_text(size=6),
        legend.position = "none")

# STOP, think
ggsave(filename = file.path(dir, 'DEdown_combinedEdgeR_woScaling.png'),
       width = 10, height = 6)
ggsave(filename = file.path(dir, 'DEdown_combinedEdgeR_Scaled.png'),
       width = 10, height = 6)

