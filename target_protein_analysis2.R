library(ggplot2)
library(ggpubr)
library(ggh4x)
library(dplyr)
library(tidyr)
library(stringr)

#############################
# Proteomics
#############################

### 1. Upload the info table from sum_orthoinfo.py
target_info <- read.csv('~/labeglo2/MS_results/390/withDBfromRNAspades/hsps_orthologes/orthologues_hsps_sum_proteinortho_with_selfblast.tsv',
                        sep = '\t') 
### 2. Upload intensities tables (after normalization) and scale it:
eve_intensities <- read.table('~/labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_eve/intensities_after_slNorm_eve.csv', 
                              header = T)
ecy_intensities <- read.table('~/labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_ecy/intensities_after_slNorm_ecy.csv', 
                              header = T)
gla_intensities <- read.table('~/labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_gla/intensities_after_slNorm_gla.csv', 
                              header = T)

data_prep <- function(data_intensities){
  row.names(data_intensities) <- data_intensities$protein_group
  data_intensities <- data_intensities[-1]
  data_intensities_scaled <- apply(data_intensities, 1, scale)
  data_intensities_scaled <- t(data_intensities_scaled)
  colnames(data_intensities_scaled) <- colnames(data_intensities)
  return(data_intensities_scaled)
}

eve_intensities <- data_prep(eve_intensities)
ecy_intensities <- data_prep(ecy_intensities)
gla_intensities <- data_prep(gla_intensities)

### 3. Upload meta file:
meta <- read.csv(file = 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/Metadata_Proteus.tsv',
                 sep = '\t')
meta$condition <- ifelse(grepl('CK1|CK2|VrK1|VrK2|LK1|LK2', meta$sample), 
                         '6C', '24C')

### 4. Connect the intensities tables and info table:
# separate the data on two parts: with the same species name and with different
true_or_false <- c()
for (row in 1:nrow(target_info)) {
  true_or_false <- c(true_or_false, grepl(target_info$species[row], 
                                          target_info$protein[row]))
}
target_info_samespecies <- target_info[true_or_false,]
target_info_diffspecies <- target_info[!true_or_false,]

# in target_info_diffspecies, leave only the proteins which can be found in transcriptome
true_or_false <- c()
for (row in 1:nrow(target_info_diffspecies)) {
  true_or_false <- c(true_or_false, 
                     grepl(target_info_diffspecies$species[row], 
                           target_info_diffspecies$protein_group[row]))
}
target_info_diff_saved <- target_info_diffspecies[true_or_false,] # will be analysed later 

connect_intens_and_info <- function(row_with_info, intens_data){
  new_df <- NULL
  pg_name <- row_with_info[2]
  intensities <- intens_data[row.names(intens_data) == as.character(pg_name),]
  for (el in 1:length(intensities)){
    if (length(intensities) == 0) {
      break
    }
    new_row <- row_with_info
    new_row[[length(new_row) + 1]] <- intensities[el]
    new_row[[length(new_row) + 1]] <- names(intensities[el])
    new_df <- rbind(new_df, new_row)
  }
  return(new_df)
}

combined_data_proteomics <- NULL
for (row in 1:nrow(target_info_samespecies)) { # data will contain only those proteins which transcripts came from the same species
  if (target_info_samespecies[row, 'species'] == 'Ecy'){
    combined_data_proteomics <- 
      rbind(combined_data_proteomics,
            connect_intens_and_info(target_info_samespecies[row,], ecy_intensities))
  }  
  else if (target_info_samespecies[row, 'species'] == 'Eve') {
    combined_data_proteomics <- 
      rbind(combined_data_proteomics,
            connect_intens_and_info(target_info_samespecies[row,], eve_intensities))
  }
  else if (target_info_samespecies[row, 'species'] == 'Gla') {
    combined_data_proteomics <- 
      rbind(combined_data_proteomics,
            connect_intens_and_info(target_info_samespecies[row,], gla_intensities))
  }
}

colnames(combined_data_proteomics)[8] <- 'intensities'
colnames(combined_data_proteomics)[9] <- 'sample'
combined_data_proteomics$condition <- ifelse(grepl('CK1|CK2|VrK1|VrK2|LK1|LK2', 
                                                   combined_data_proteomics$sample), 
                                             '6C', '24C')

###########################################
# Transcriptomics
###########################################

### 1. Choose the species:
species <- 'Ecy'
species <- 'Eve'
species <- 'Gla'

### 2. Upload counts table from quantification of transcriptomic data
dir <- paste0('~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/', 
              species, '/counts') # specify path to your samples

countFiles <- list.files(dir, full.names = T)
countFiles
counts <- lapply(countFiles, function(countsFile) {
  f <- file.path(countsFile, 'quant.sf')
  read.table(f, sep="\t", header=1, row.names = 1, 
             stringsAsFactors = F, comment.char = "#")
})
### 3. Take TPMs (before I took NumReads, but here is better to take normalized data)
counts <- lapply(counts, function(countsTable) countsTable[, 3, drop=F])

counts <- do.call(cbind, counts)
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3', '24C_rep4',
                  '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4') # Ecy, Eve
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3',
                  '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4') # Gla
colnames(counts) <- sample_names

# Need to be write properly (not to repeat every time)
counts_ecy <- counts
counts_eve <- counts
counts_gla <- counts

### 4. Prepare counts and scale the data:
data_prep_counts <- function(data_counts){
  data_counts_scaled <- apply(data_counts, 1, scale)
  data_counts_scaled <- t(data_counts_scaled)
  colnames(data_counts_scaled) <- colnames(data_counts)
  return(data_counts_scaled)
}

counts_ecy <- data_prep_counts(counts_ecy)
counts_eve <- data_prep_counts(counts_eve)
counts_gla <- data_prep_counts(counts_gla)

### 5. Connect counts and the info table:

connect_counts_and_info <- function(row_with_info, counts_data, species_name){
  new_df <- NULL
  #  if (grepl(species_name, row_with_info[1]) == TRUE) {
  #contig_name <- sub(paste0('_', species_name), '', row_with_info[1])
  #contig_name <- sub('.p[1-9]*', '', contig_name)
  contig_name <- sub('.p[1-9]*_...$', '', row_with_info[1])
  counts_list <- 
    counts_data[row.names(counts_data) == as.character(contig_name),]
  for (el in 1:length(counts_list)){
    if (length(counts_list) == 0) {
      break
    }
    new_row <- row_with_info
    new_row[[length(new_row) + 1]] <- counts_list[el]
    new_row[[length(new_row) + 1]] <- names(counts_list[el])
    new_df <- rbind(new_df, new_row)
  }
  #  }
  return(new_df)
}

combined_data_transc <- NULL # data will contain only those proteins which transcripts came from the same species
for (row in 1:nrow(target_info_samespecies)) {
  if (target_info_samespecies[row, 'species'] == 'Ecy'){
    combined_data_transc <- 
      rbind(combined_data_transc,
            connect_counts_and_info(target_info_samespecies[row,], counts_ecy, 'Ecy'))
  }  
  else if (target_info_samespecies[row, 'species'] == 'Eve') {
    combined_data_transc <- 
      rbind(combined_data_transc,
            connect_counts_and_info(target_info_samespecies[row,], counts_eve, 'Eve'))
  }
  else if (target_info_samespecies[row, 'species'] == 'Gla') {
    combined_data_transc <- 
      rbind(combined_data_transc,
            connect_counts_and_info(target_info_samespecies[row,], counts_gla, 'Gla'))
  }
}

colnames(combined_data_transc)[8] <- 'counts'
colnames(combined_data_transc)[9] <- 'sample'
combined_data_transc$condition <- ifelse(grepl('6C', combined_data_transc$sample), 
                                         '6C', '24C')

#################################
# Proteomics + Transcriptomics
#################################
colnames(combined_data_proteomics) <- c(colnames(combined_data_proteomics)[1:7],
                                        "values", 
                                        colnames(combined_data_proteomics)[9:10])
combined_data_proteomics$method <- 'MS/MS'
colnames(combined_data_transc) <- c(colnames(combined_data_transc)[1:7],
                                    "values", 
                                    colnames(combined_data_transc)[9:10])
combined_data_transc$method <- 'RNAseq'
#combined_data <- rbind(combined_data_prot_samespecies, combined_data_transc)
combined_data_proteomics$unique_name <- sprintf('%s_%s', 
                                                combined_data_proteomics$protein,
                                                combined_data_proteomics$species)
combined_data_transc$unique_name <- sprintf('%s_%s', 
                                            combined_data_transc$protein,
                                            combined_data_transc$species)
combined_data <- NULL
for (row in 1:nrow(combined_data_proteomics)) {
  if (any(grepl(combined_data_proteomics[row,][12], 
                combined_data_transc$unique_name))) {
    rows2add <- combined_data_transc[combined_data_transc$unique_name == as.character(combined_data_proteomics[row,][12]),]
    combined_data <- rbind(combined_data, combined_data_proteomics[row,],
                           rows2add)
  }
  else {
    print(combined_data_proteomics[row,][12])
  }
}
combined_data <- unique(combined_data)
# 69 proteins out of 71 left: NODE_16065_length_1920_cov_60.704436_g10500_i0.p1_Eve_Eve
# and NODE_6760_length_3076_cov_1784.947803_g3984_i0.p1_Gla_Gla seems to be not here
# They were filtered out (perhaps, all zeros in the raw for intensities)

#####################################################################################
### Let's work with proteins which have the matched species name inside protein group
#####################################################################################
target_info_diff_saved
target_info_diff_saved$pg <- target_info_diff_saved$protein_group
target_info_diff_saved_sep <- separate(target_info_diff_saved, protein_group, 
                                       ';', into=paste0('key', 1:30), 
                         fill='right')
target_info_diff_saved_long <- pivot_longer(target_info_diff_saved_sep, 
                                            starts_with('key'),
                                            names_to = NULL, 
                                            values_to = 'protein_from_pg') %>%
                                            filter(!is.na(protein_from_pg))
true_or_false <- c()
for (row in 1:nrow(target_info_diff_saved_long)) {
  true_or_false <- c(true_or_false, grepl(target_info_diff_saved_long$species[row], 
                                          target_info_diff_saved_long$protein_from_pg[row]))
}
target_diff_saved_filt <- target_info_diff_saved_long[true_or_false,]

td1 <- subset(target_diff_saved_filt, 
              protein != 'NODE_3822_length_4220_cov_2696.563654_g2376_i0.p1_Eve')

td2 <- subset(target_diff_saved_filt, 
              protein_from_pg == 'NODE_5378_length_3701_cov_2676.348850_g3383_i0.p1_Ecy')

target_diff_saved_filt <- rbind(td1, td2)
target_diff_saved_filt$protein <- target_diff_saved_filt$protein_from_pg
target_diff_saved_filt$protein_group <- target_diff_saved_filt$pg
target_diff_saved_filt_test <- target_diff_saved_filt[,c(1, 9, 2:6)]
########################################################################################
## Connect obtained target data (target_diff_saved_filt) with intensities and counts 
combined_data_prot_2 <- NULL
for (row in 1:nrow(target_diff_saved_filt_test)) { # data will contain only those proteins which transcripts came from the same species
  if (target_diff_saved_filt_test[row, 'species'] == 'Ecy'){
    combined_data_prot_2 <- 
      rbind(combined_data_prot_2,
            connect_intens_and_info(target_diff_saved_filt_test[row,], ecy_intensities))
  }  
  else if (target_diff_saved_filt_test[row, 'species'] == 'Eve') {
    combined_data_prot_2 <- 
      rbind(combined_data_prot_2,
            connect_intens_and_info(target_diff_saved_filt_test[row,], eve_intensities))
  }
  else if (target_diff_saved_filt_test[row, 'species'] == 'Gla') {
    combined_data_prot_2 <- 
      rbind(combined_data_prot_2,
            connect_intens_and_info(target_diff_saved_filt_test[row,], gla_intensities))
  }
}

colnames(combined_data_prot_2)[8] <- 'intensities'
colnames(combined_data_prot_2)[9] <- 'sample'
combined_data_prot_2$condition <- ifelse(grepl('CK1|CK2|VrK1|VrK2|LK1|LK2', 
                                        combined_data_prot_2$sample), 
                                             '6C', '24C')

combined_data_transc_2 <- NULL # data will contain only those proteins which transcripts came from the same species
for (row in 1:nrow(target_diff_saved_filt_test)) {
  if (target_diff_saved_filt_test[row, 'species'] == 'Ecy'){
    combined_data_transc_2 <- 
      rbind(combined_data_transc_2,
            connect_counts_and_info(target_diff_saved_filt_test[row,], counts_ecy, 'Ecy'))
  }  
  else if (target_diff_saved_filt_test[row, 'species'] == 'Eve') {
    combined_data_transc_2 <- 
      rbind(combined_data_transc_2,
            connect_counts_and_info(target_diff_saved_filt_test[row,], counts_eve, 'Eve'))
  }
  else if (target_diff_saved_filt_test[row, 'species'] == 'Gla') {
    combined_data_transc_2 <- 
      rbind(combined_data_transc_2,
            connect_counts_and_info(target_diff_saved_filt_test[row,], counts_gla, 'Gla'))
  }
}

colnames(combined_data_transc_2)[8] <- 'counts'
colnames(combined_data_transc_2)[9] <- 'sample'
combined_data_transc_2$condition <- ifelse(grepl('6C', combined_data_transc_2$sample), 
                                         '6C', '24C')
### Combine intensity data and counts
colnames(combined_data_prot_2) <- c(colnames(combined_data_prot_2)[1:7],
                                        "values", 
                                        colnames(combined_data_prot_2)[9:10])
combined_data_prot_2$method <- 'MS/MS'
colnames(combined_data_transc_2) <- c(colnames(combined_data_transc_2)[1:7],
                                    "values", 
                                    colnames(combined_data_transc_2)[9:10])
combined_data_transc_2$method <- 'RNAseq'
#combined_data <- rbind(combined_data_prot_samespecies, combined_data_transc)
combined_data_prot_2$unique_name <- sprintf('%s_%s', 
                                            combined_data_prot_2$protein,
                                            combined_data_prot_2$species)
combined_data_transc_2$unique_name <- sprintf('%s_%s', 
                                              combined_data_transc_2$protein,
                                              combined_data_transc_2$species)
combined_data_2 <- NULL
for (row in 1:nrow(combined_data_prot_2)) {
  if (any(grepl(combined_data_prot_2[row,][12], 
                combined_data_transc_2$unique_name))) {
    rows2add <- combined_data_transc_2[combined_data_transc_2$unique_name == as.character(combined_data_prot_2[row,][12]),]
    combined_data_2 <- rbind(combined_data_2, combined_data_prot_2[row,],
                           rows2add)
  }
  else {
    print(combined_data_prot_2[row,][12])
  }
}
combined_data_2 <- unique(combined_data_2)
##################################################
#### Combine combine_data and combine_data_2

combined_data_full <- rbind(combined_data, combined_data_2)

### To add information about significance of the difference
total_transc <- NULL
total_proteins <- NULL
all_species <- c('Eve', 'Ecy', 'Gla')
for (species in 1:length(all_species)) {
  sp <- all_species[species]
  transcr <- 
    read.csv(paste0('~/labeglo2/proteome_transcr_comparision/', 
                    sp, '_transcr_24vs6_all.csv'), sep = '\t')
  transcr$unique_name <- paste0(transcr$contig, '_', sp)
  total_transc <- rbind(total_transc, transcr)
  proteins <- 
    read.csv(paste0('~/labeglo2/proteome_transcr_comparision/', 
                    sp, 
                    '_AllProteins_24vs6after_proteinGroups_separatelyAnalyzed.csv'),
             sep = '\t', header = T)
  proteins$unique_name <- paste0(proteins$protein, '_', sp)
  total_proteins <- rbind(total_proteins, proteins)
}

combined_data_full$pg_unique_name <- paste0(combined_data_full$protein_group, '_', 
                                            combined_data_full$species)
combined_data_full$contig_unique_name <- sub('\\.p[1-9]*_...$', '', 
                                             combined_data_full$protein)
combined_data_full$contig_unique_name <- paste0(combined_data_full$contig_unique_name, '_', 
                                           combined_data_full$species)

combined_data_full <- cbind(combined_data_full,
                       total_proteins[match(combined_data_full$pg_unique_name, total_proteins$unique_name),][6])
combined_data_full <- cbind(combined_data_full, 
                       total_transc[match(combined_data_full$contig_unique_name, total_transc$unique_name),][6])

#total_proteins[total_proteins$protein == 'NODE_9313_length_2618_cov_250.297392_g5535_i0.p1_Gla',]

combined_data_full$sign2plot <- 
  ifelse(combined_data_full$FDR_recalc < 0.05 & combined_data_full$padj < 0.05,
         '< 0.05', 
         ifelse(combined_data_full$padj > 0.05 & combined_data_full$method == 'RNAseq',
                '> 0.05', 
                ifelse(combined_data_full$FDR_recalc > 0.05 & combined_data_full$method == 'MS/MS',
                       '> 0.05', '< 0.05')))
combined_data_full$method <- factor(combined_data_full$method, 
                               levels = c('RNAseq', 'MS/MS'))

## here is a handling with special case (NODE_9313_length_2618_cov_250.297392_g5535_i0.p1_Gla)
# we have values from both proteomics and transcriptomics but I do not have the info about 
# p-value from proteomics. But obviously there is no difference -> put <0.05:
combined_data_full_gla_ecy70const <- subset(combined_data_full, 
                                            protein == 'NODE_9313_length_2618_cov_250.297392_g5535_i0.p1_Gla')
combined_data_full_gla_ecy70const$sign2plot <- 
  ifelse(combined_data_full_gla_ecy70const$method == 'RNAseq',
         '< 0.05', '> 0.05')
combined_data_full <- subset(combined_data_full, 
                             protein != 'NODE_9313_length_2618_cov_250.297392_g5535_i0.p1_Gla')
combined_data_full <- rbind(combined_data_full, combined_data_full_gla_ecy70const)
###########################
## To plot
###########################
to_plot <- subset(combined_data_full, orthogroup == 8)
to_plot$annotation <- sub(' isoform[^,]*,?', '', to_plot$annotation)
to_plot$annotation <- sub(' partial', '', to_plot$annotation)
to_plot$annotation <- sub('PREDICTED: ', '', to_plot$annotation)
to_plot$annotation <- sub('-like', '', to_plot$annotation)
var_width <- 20
to_plot <- mutate(to_plot, 
                  pretty_varname = str_wrap(to_plot$annotation, 
                                            width = var_width))
to_plot$all_labels <- sprintf("%s%s|%s, %s",
                              to_plot$species, to_plot$protein_group,
                              to_plot$pretty_varname, to_plot$hsp70_type)
#to_plot$all_labels <- sprintf("%s%s|%s, %s",
#                              to_plot$species, to_plot$protein_group,
#                              to_plot$annotation, to_plot$hsp70_type)
f <- function(x) {
  sapply(strsplit(x, '|', fixed=T), `[`, 2)
}
to_plot$species <- factor(to_plot$species, levels = c('Eve', 'Gla', 'Ecy'), 
                          labels = c('E.verrucosus', 'G.lacustris', 
                                     'E.cyaneus'))
to_plot$species_italic <- sprintf('italic(%s)', to_plot$species)

#to_plot <- na.omit(to_plot)
to_plot$condition <- factor(to_plot$condition, levels = c('24C', '6C'),
                            labels = c('24 °C', '6 °C'))
to_plot$species_italic <- factor(to_plot$species_italic,
                                 levels = c('italic(E.cyaneus)',
                                            'italic(E.verrucosus)',
                                            'italic(G.lacustris)'))
to_plot$species_italic <- factor(to_plot$species_italic,
                                 levels = c('italic(E.verrucosus)',
                                            'italic(G.lacustris)',
                                            'italic(E.cyaneus)'))
to_plot$species_italic <- factor(to_plot$species_italic,
                                 levels = c('italic(E.cyaneus)',
                                            'italic(G.lacustris)',
                                            'italic(E.verrucosus)'))
##### if do not want ', *' in the end of annotations:
to_plot$all_labels <- sub(', \\*', '', to_plot$all_labels)
#####

ggplot(to_plot, aes(x = method, y = values, color = condition)) +
  geom_point(position=position_jitterdodge(dodge.width=.8)) +
  geom_boxplot(aes(fill = condition, linetype = sign2plot),
              outlier.alpha = 0, alpha = 0.4) +
             # ,
              #linetype = 'dashed') 
  facet_nested_wrap(species_italic ~ all_labels, 
                    labeller = labeller(all_labels = as_labeller(f), 
                                        species_italic = label_parsed), 
                    scales = 'free', ncol = 5) +
  scale_color_manual('Condition:', values = c('#ca0020', '#0571b0')) +
  scale_fill_manual('Condition:', values = c('#ca0020', '#0571b0')) +
  scale_linetype('adj. p-value:') +
  xlab('') +
  ylab('Scaled absolute values') +
  theme_bw() 
#  theme(strip.text = element_text(size=5.5))

dir_to_save <- '/home/polina/labeglo2/MS_results/390/withDBfromRNAspades/hsps_orthologes'
ggsave(file.path(dir_to_save, 'og32_proteinortho_proteinsWITHtranscripts.png'),
       #scale = 1.2) 
       width = 3, height = 3)
ggsave(file.path(dir_to_save, 'og32_proteinortho_proteinsWITHtranscripts.pdf'),
       #scale = 1.2)
       width = 3.3, height = 3)

# og1 - width = 11
# og2 - width = 11, height = 3
# og21 - width = 5, height = 3 (two tiles)
# og25 - width = 3, height = 3 (one tile)






