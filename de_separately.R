# DE analysis for the each group with missing values (without NAs, with 3 NAs 
# (3 TMT batch), 8 NAs (1 TMT batch), 8 NAs (2 TMT batch), 11 NAs (1 and 3 TMT 
# batches), 11 NAs (1 and 2 TMT batches))

library(tidyverse) 
library(limma) 
library(edgeR)
library(ggfortify) 

species <- 'Ecy'

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

meta <- meta_upload(path2meta, species)

# if you want to take all 6C controls as one control group:
meta$condition <- ifelse(grepl('CK1|CK2|VrK1|VrK2|LK1|LK2', meta$sample), 
                         '6C', '24C')

### 2. Data uploading (proteinGroups file from MaxQuant)
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_eve/'
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_ecy/'
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_gla/'
proteinGroups_file <- 'proteinGroups_wo_cont_more2pept.txt' # take file with proteinGroups with 2 or more peptides quantified
dat_init <- read.csv(file.path(dir, proteinGroups_file), sep = '\t', header = T, 
                     check.names = F) 
pep_annot <- read.csv(file.path(dir, 'annot_protein_groups_ecy.csv'), sep = '\t',
                      header = T)

select_data <- function(meta_data, proteinGroups_data){
  dat_ <- proteinGroups_data[,c('Protein IDs', meta_data$measure)]
  rownames(dat_) <- dat_$`Protein IDs`
  dat_ <- dat_[-1]
  colnames(dat_) <- meta_data$sample
  dat_[dat_ == 0] <- NA
  return(dat_)
}
dat <- select_data(meta, dat_init)

### 3. Get rid of the outliers and samples, which you do not want to analyze
dat <- dat[,!grepl('VrK1_390_3|VrK1_390_4|pool', colnames(dat))]
meta <- subset(meta, sample != 'VrK1_390_3' & sample != 'VrK1_390_4' &
                 !grepl('pool', sample))

#'L1_390_3|L1_390_4'

### 4. Remove rows with all NAs in at least one of the condition

remove_completeNAs_inCond <- function(data2clean, metafile){
  cond <- as.factor(metafile$condition)
  ixs <- order(cond)
  data2clean <- data2clean[, ixs]
  cond <- cond[ixs]
  
  for (level in levels(cond)) {
    col_ixs <- cond == level
    data2clean <- data2clean[rowSums(is.na(data2clean[,col_ixs])) < sum(col_ixs),]
  }
  return(data2clean)
}
data_raw <- remove_completeNAs_inCond(dat, meta)

### 5. Visualize the raw data
boxplot(log2(data_raw), col = c(rep('blue', 8), rep('red', 5), rep('green', 6)), 
        notch = TRUE, main = 'RAW data: 24C (blue), 6C_after (red), 6C_before (green)',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_raw), col = c(rep('blue', 8), rep('red', 5), rep('green', 6)), 
              main = 'Raw data', legend=F) # Eve, resolution -> 800/553

boxplot(log2(data_raw), col = c(rep('blue', 9), rep('red', 8)), 
        notch = TRUE, main = 'RAW data: 24C (blue), 6C_after (red), 6C_before (green)',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_raw), col = c(rep('blue', 9), rep('red', 8)), 
              main = 'Raw data', legend=F) # Ecy

boxplot(log2(data_raw), col = c(rep('blue', 9), rep('red', 9), rep('green', 4)), 
        notch = TRUE, main = 'RAW data: 24C (blue), 6C_after (red), 6C_before (green)',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_raw), col = c(rep('blue', 9), rep('red', 9), 
                                      rep('green', 4)), 
              main = 'Raw data', legend=F) # Gla

### 6. Sample loading normalization

# Eve
exp1_raw <- data_raw[,subset(meta, experiment == 6)$sample]
exp2_raw <- data_raw[,subset(meta, experiment == 7)$sample]
exp3_raw <- data_raw[,subset(meta, experiment == 8)$sample]

# Ecy
exp1_raw <- data_raw[,subset(meta, experiment == 4)$sample]
exp2_raw <- data_raw[,subset(meta, experiment == 5)$sample]

# Gla
exp1_raw <- data_raw[,subset(meta, experiment == 1)$sample]
exp2_raw <- data_raw[,subset(meta, experiment == 2)$sample]
exp3_raw <- data_raw[,subset(meta, experiment == 3)$sample]

target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T), 
               colSums(exp3_raw, na.rm = T)),
               na.rm = T)
target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T)),
               na.rm = T)
norm_facs <- target / colSums(exp1_raw, na.rm = T)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw, na.rm = T)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3_raw, na.rm = T)
exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")

data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl)
data_sl <- cbind(exp1_sl, exp2_sl)

boxplot(log2(data_sl), col = c(rep(c('red', 'blue'), each = 8), rep('green', 3)), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nExp1 (red), Exp2 (blue), Exp3 (green)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_sl), col = c(rep(c('red', 'blue'), each = 8), rep('green', 3)),
              main = "SL normalization", legend = F) # Eve, resolution -> 800/553

boxplot(log2(data_sl), col = c(rep('red', 9), rep('green', 8)), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nExp1 (red), Exp2 (blue), Exp3 (green)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity') # Ecy 
plotDensities(log2(data_sl), col = c(rep('red', 9), rep('green', 8)),
              main = "SL normalization", legend = F) # Ecy

boxplot(log2(data_sl), col = c(rep(c('red', 'blue'), each = 9), rep('green', 4)), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nExp1 (red), Exp2 (blue), Exp3 (green)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity') # Gla
plotDensities(log2(data_sl), col = c(rep(c('red', 'blue'), each = 9), 
                                     rep('green', 4)),
              main = "SL normalization", legend = F) # Gla

### 7. Vizualize normalized data
## 7.1 MDS-plot with groups colored by TMT batch
plotMDS(log2(data_sl), col = c(rep(c('red', 'blue'), each = 8), rep('green', 3)), 
        main = "SL clusters grouped by TMT experiment") # Eve
plotMDS(log2(data_sl), col = c(rep('red', 9), rep('green', 8)), 
        main = "SL clusters grouped by TMT experiment") # Ecy
plotMDS(log2(data_sl), col = c(rep(c('red', 'blue'), each = 9), rep('green', 4)), 
        main = "SL clusters grouped by TMT experiment") # Gla

## 7.2 MDS-plot with groups colored by condition:
png(file.path(dir, 'Gla_MDS_wIMBR_woImput_woPools.png'),
    width = 1200, height = 1000, pointsize = 28)
col_vec <- ifelse(grepl('pool', colnames(data_sl)), 'green',
                  ifelse(grepl('VrK1|CK1|LK1', colnames(data_sl)), 'orange',
                         ifelse(grepl('Vr1|C1|L1', colnames(data_sl)), 'blue', 'red')))

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
data_sl_mult$protein_group <- row.names(data_sl_mult)
data_sl_mult_ <- data.frame(protein_group = data_sl_mult$protein_group, 
                            data_sl_mult[1:length(data_sl_mult)-1])
write.table(data_sl_mult_, file.path(dir, 'intensities_after_slNorm_ecy.csv'), 
            sep = '\t', quote = F, row.names = F)

#########
# EdgeR 
#########
### 9. Separate the data on the subset with the same NA pattern and perform EdgeR analysis

# here are two functions to perform de analysis
edger_analysis <- function(datatode, group, time_group, cond_vector, 
                           protein_annotation){
  
  edge_object <- DGEList(counts = datatode, group = group)
  if (length(unique(time_group)) > 1){
    design <- model.matrix(~0 + group + time_group)
  }
  else {
    design <- model.matrix(~0 + group)
  }
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

edger_analysis_wrap <- function(norm_data, protein_annotation, 
                                cond2compare, metafile){
  cond2compare_ <- paste(cond2compare, collapse = '|') # to make it readable for grepl
  samples2take <- metafile[grepl(cond2compare_, metafile$condition),]$sample
  samples2take <- paste(samples2take, collapse = '|')
  data2DE <- norm_data[,grepl(samples2take, colnames(norm_data))] 
  
  sample_order <- match(colnames(data2DE), metafile$sample)
  group <- metafile$condition[sample_order]
  
  #group <- ifelse(grepl('VrK1|VrK2|CK1|CK2|LK1|LK2', 
  #                      colnames(norm_data)), '6C', '24C')
  group_time <- as.factor(ifelse(grepl('VrK1|CK1|LK1', colnames(data2DE)), 25, 26))
  
  print("Samples to be analyzed:")
  print(colnames(data2DE))
  
  return(edger_analysis(data2DE, group, group_time, cond2compare, 
                        protein_annotation))
}

# wo NAs
data_wo_na <- na.omit(data_irs)
SignProteinInfo_0 <- edger_analysis_wrap(norm_data = data_wo_na, 
                                         protein_annotation = pep_annot, 
                                         cond2compare = c("6C", 
                                                          "24C"), 
                                         metafile = meta)

# 3 NAs/row - Eve
data_3 <- data_irs[rowSums(is.na(data_irs)) == 3,]
data_3 <- data_3[,colSums(is.na(data_3))<nrow(data_3)]
SignProteinInfo_3 <- edger_analysis_wrap(norm_data = data_3, 
                                         protein_annotation = pep_annot, 
                                         cond2compare = c("6C", 
                                                          "24C"), 
                                         metafile = meta)

# 8 NAs/row (the 1st + 3rd TMT batches) - Eve
data_8_1 <- data_irs[rowSums(is.na(data_irs)) == 8,]
data_8_1 <- data_8_1[,grepl('_1$|_2$|_3$|Vr1_390_7|Vr1_390_8|VrK1_390_8', 
                          colnames(data_8_1))]
data_8_1 <- na.omit(data_8_1)
SignProteinInfo_8_1 <- edger_analysis_wrap(norm_data = data_8_1, 
                                           protein_annotation = pep_annot, 
                                           cond2compare = c("6C", 
                                                            "24C"), 
                                           metafile = meta)

# 8 NAs/row (the 2nd + 3rd TMT batches)
data_8_2 <- data_irs[rowSums(is.na(data_irs)) == 8,]
data_8_2 <- data_8_2[,grepl('_4$|_5$|_6$|VrK1_390_7|Vr1_390_7|Vr1_390_8|VrK1_390_8', 
                            colnames(data_8_2))]
data_8_2 <- na.omit(data_8_2)
SignProteinInfo_8_2 <- edger_analysis_wrap(norm_data = data_8_2, 
                                           protein_annotation = pep_annot, 
                                           cond2compare = c("6C", 
                                                            "24C"), 
                                           metafile = meta)

# 11 NAs/row (the 1st TMT batch)
data_11_1 <- data_irs[rowSums(is.na(data_irs)) == 11,]
data_11_1 <- data_11_1[,!grepl('Vr1_390_7|Vr1_390_8|VrK1_390_8', 
                             colnames(data_11_1))]
data_11_1 <- data_11_1[,grepl('_1$|_2$|_3$', 
                            colnames(data_11_1))]
data_11_1 <- na.omit(data_11_1)
SignProteinInfo_11_1 <- edger_analysis_wrap(norm_data = data_11_1, 
                                            protein_annotation = pep_annot, 
                                            cond2compare = c("6C", 
                                                             "24C"), 
                                            metafile = meta)

# 11 NAs/row (the 2nd TMT batch)
data_11_2 <- data_irs[rowSums(is.na(data_irs)) == 11,]
data_11_2 <- data_11_2[,!grepl('Vr1_390_7|Vr1_390_8|VrK1_390_8', 
                             colnames(data_11_2))]
data_11_2 <- data_11_2[,grepl('_4$|_5$|_6$|VrK1_390_7', 
                            colnames(data_11_2))]
data_11_2 <- na.omit(data_11_2)
SignProteinInfo_11_2 <- edger_analysis_wrap(norm_data = data_11_2, 
                                            protein_annotation = pep_annot, 
                                            cond2compare = c("6C", 
                                                             "24C"), 
                                            metafile = meta)

# 9 NAs/row (the 2nd TMT batch) - Ecy
data_9 <- data_irs[rowSums(is.na(data_irs)) == 9,]
data_9 <- data_9[,meta[meta$experiment == 5,]$sample]
SignProteinInfo_9 <- edger_analysis_wrap(norm_data = data_9, 
                                         protein_annotation = pep_annot, 
                                         cond2compare = c("6C", 
                                                          "24C"), 
                                         metafile = meta)

# 8 NAs/row (the 1st TMT batch) - Ecy
data_8 <- data_irs[rowSums(is.na(data_irs)) == 8,]
data_8 <- data_8[,meta[meta$experiment == 4,]$sample]
SignProteinInfo_8 <- edger_analysis_wrap(norm_data = data_8, 
                                         protein_annotation = pep_annot, 
                                         cond2compare = c("6C", 
                                                          "24C"), 
                                         metafile = meta)

# 4 NAs/row - Gla
data_4 <- data_irs[rowSums(is.na(data_irs)) == 4,]
data_4 <- data_4[,meta[meta$experiment != 3,]$sample]
SignProteinInfo_4 <- edger_analysis_wrap(norm_data = data_4, 
                                         protein_annotation = pep_annot, 
                                         cond2compare = c("6C", 
                                                          "24C"), 
                                         metafile = meta)

# 9 NAs/row (the 1st TMT batch) - Gla
data_9_1 <- data_irs[rowSums(is.na(data_irs)) == 9,]
data_9_1 <- data_9_1[,meta[meta$experiment != 2,]$sample]
data_9_1 <- na.omit(data_9_1)
SignProteinInfo_9_1 <- edger_analysis_wrap(norm_data = data_9_1, 
                                           protein_annotation = pep_annot, 
                                           cond2compare = c("6C", 
                                                            "24C"), 
                                           metafile = meta)

# 9 NAs/row (the 2nd TMT batch) - Gla
data_9_2 <- data_irs[rowSums(is.na(data_irs)) == 9,]
data_9_2 <- data_9_2[,meta[meta$experiment != 1,]$sample]
data_9_2 <- na.omit(data_9_2)
SignProteinInfo_9_2 <- edger_analysis_wrap(norm_data = data_9_2, 
                                           protein_annotation = pep_annot, 
                                           cond2compare = c("6C", 
                                                            "24C"), 
                                           metafile = meta)

# 9 NAs/row (the 2nd TMT batch) - Gla
data_13_1 <- data_irs[rowSums(is.na(data_irs)) == 13,]
data_13_1 <- data_13_1[,meta[meta$experiment != 3 & meta$experiment != 2,]$sample]
data_13_1 <- na.omit(data_13_1)
SignProteinInfo_13_1 <- edger_analysis_wrap(norm_data = data_13_1, 
                                            protein_annotation = pep_annot, 
                                            cond2compare = c("6C", 
                                                             "24C"), 
                                            metafile = meta)

# 9 NAs/row (the 1st TMT batch) - Gla
data_13_2 <- data_irs[rowSums(is.na(data_irs)) == 13,]
data_13_2 <- data_13_2[,meta[meta$experiment != 3 & meta$experiment != 1,]$sample]
data_13_2 <- na.omit(data_13_2)
SignProteinInfo_13_2 <- edger_analysis_wrap(norm_data = data_13_2, 
                                            protein_annotation = pep_annot, 
                                            cond2compare = c("6C", 
                                                             "24C"), 
                                            metafile = meta)

### 10. Combine all subsets in one and calculate the adjusted p-values
SPI_all <- rbind(SignProteinInfo_0, SignProteinInfo_3, SignProteinInfo_8_1,
                 SignProteinInfo_8_2, SignProteinInfo_11_1, SignProteinInfo_11_2)

SPI_all <- rbind(SignProteinInfo_0, SignProteinInfo_8,
                 SignProteinInfo_9) # Ecy

SPI_all <- rbind(SignProteinInfo_0, SignProteinInfo_4,
                 SignProteinInfo_9_1, SignProteinInfo_9_2,
                 SignProteinInfo_13_1, SignProteinInfo_13_2) # Gla

SPI_all$FDR_recalc <- p.adjust(SPI_all$pvalue, method = "fdr")

SPI_all_sign <- subset(SPI_all, FDR_recalc < 0.05)

nrow(subset(SPI_all_sign, logFC > 0)) # UP
nrow(subset(SPI_all_sign, logFC < 0)) # DOWN

write.table(SPI_all, 
            file = file.path(paste0('~/labeglo2/proteome_transcr_comparision/',
                                    species, '_AllProteins_24vs6after_proteinGroups_separatelyAnalyzed.csv')),
            sep = '\t', quote = F, row.names = F)

########################################
# 11. Draw boxplots for the DE proteins
########################################
to_plot <- SPI_all_sign
#to_plot <- subset(SignProteinInfo_0, FDR < 0.05)
res_sign_fdr_up <- subset(to_plot, logFC > 0)
res_sign_fdr_down <- subset(to_plot, logFC < 0)
dat_mult_stand_fdr <- apply(data_sl_mult, 1, scale) # if take multiplied and scaled data to visualize
dat_mult_stand_fdr <- t(dat_mult_stand_fdr) # --||--
#dat_mult_stand_fdr <- data_sl # if take "naked" data after normalization
colnames(dat_mult_stand_fdr) <- colnames(data_sl_mult)
dat_sign_fdr <- dat_mult_stand_fdr[row.names(dat_mult_stand_fdr) %in% 
                                     to_plot$protein,]

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
dat_sign_fdr_up_long$Condition <- ifelse(grepl('Vr1|C1|L1', dat_sign_fdr_up_long$name), 
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
             ncol = 6) +
  ylab('Scaled intensities') +
  xlab('') +
  theme_bw() +
  theme(strip.text = element_text(size=8),
        legend.position = "none")

# STOP, think
ggsave(filename = file.path(dir, 'DEup_combinedEdgeR_Scaled_24vs6C.png'),
       width = 9.6, height = 6)
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
dat_sign_fdr_down_long$Condition <- ifelse(grepl('Vr1|C1|L1', 
                                                 dat_sign_fdr_down_long$name), 
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
  ylab('Scaled intensities') +
  xlab('') +
  theme_bw() +
  theme(strip.text = element_text(size=8),
        legend.position = "none")

# STOP, think
ggsave(filename = file.path(dir, 'DEdown_combinedEdgeR_Scaled_24vs6C.png'),
       width = 10.8, height = 6)
ggsave(filename = file.path(dir, 'DEdown_combinedEdgeR_woScaling.png'),
       width = 10, height = 6)


