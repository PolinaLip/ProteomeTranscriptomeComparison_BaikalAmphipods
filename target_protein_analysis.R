library(ggplot2)
library(ggpubr)
library(ggh4x)

#############################
# Proteomics
#############################

### 1. Upload the info table from sum_orthoinfo.py
target_info <- read.csv('~/labeglo2/MS_results/390/withDBfromRNAspades/hsps_orthologes/orthologues_hsps_sum_updated.tsv',
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
for (row in 1:nrow(target_info)) {
  if (target_info[row, 'species'] == 'Ecy'){
    combined_data_proteomics <- 
      rbind(combined_data_proteomics,
        connect_intens_and_info(target_info[row,], ecy_intensities))
  }  
  else if (target_info[row, 'species'] == 'Eve') {
    combined_data_proteomics <- 
      rbind(combined_data_proteomics,
            connect_intens_and_info(target_info[row,], eve_intensities))
  }
  else if (target_info[row, 'species'] == 'Gla') {
    combined_data_proteomics <- 
      rbind(combined_data_proteomics,
            connect_intens_and_info(target_info[row,], gla_intensities))
  }
}

colnames(combined_data_proteomics)[8] <- 'intensities'
colnames(combined_data_proteomics)[9] <- 'sample'
combined_data_proteomics$condition <- ifelse(grepl('CK1|CK2|VrK1|VrK2|LK1|LK2', 
                                             combined_data_proteomics$sample), 
                                             '6C', '24C')

### 5. Plot the intensities of the found proteins:
to_plot <- subset(combined_data_proteomics, orthogroup == 'OG0000002')
to_plot$annotation <- sub(' isoform[^,]*,?', '', to_plot$annotation)
to_plot$annotation <- sub(' partial', '', to_plot$annotation)
to_plot$annotation <- sub('PREDICTED: ', '', to_plot$annotation)
to_plot$annotation <- sub('-like', '', to_plot$annotation)
to_plot$all_labels <- sprintf("%s %s|%s, %s", to_plot$species, to_plot$protein_group,
                              to_plot$annotation, to_plot$hsp70_type)

f <- function(x) {
  sapply(strsplit(x, '|', fixed=T), `[`, 2)
}
my_labeles <- c('Eve' = 'E. verrucosus', 
                'Ecy' = 'E. cyaneus', 
                'Gla' = 'G. lacustris')
my_colors <- c('#7570b3', '#1b9e77', '#d95f02')
ggplot(to_plot, aes(x = condition, y = intensities)) +
  geom_boxplot(aes(color = species, fill = species), 
               outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color = species)) +
  scale_color_manual('Species', values = my_colors,
                     labels = my_labeles) +
  scale_fill_manual('Species', values = my_colors,
                    labels = my_labeles) +
  facet_wrap(~all_labels, labeller = as_labeller(f), ncol = 4) +
  ylab('Scaled intensities') +
  xlab('') +
  theme_bw() +
  theme(legend.text = element_text(face = 'italic', size = 11),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size=7))

ggsave('~/labeglo2/MS_results/390/withDBfromRNAspades/hsps_orthologes/og19.png',
       #scale = 1.2) 
       width = 6, height = 2.7)
ggsave('~/labeglo2/MS_results/390/withDBfromRNAspades/hsps_orthologes/og19.pdf',
       #scale = 1.2)
       width = 6, height = 2.7)
# og0 - width = 11, height = 7
# og2 - default, scale = 1.2
# og1 - width = 11, height = 6
# og5 - width = 10, height = 2.7
# og7 - width = 9, height = 2.5, strip.text = element_text(size=6)
# og11, og13 -  width = 9, height = 2.7
# og19 - width = 6, height = 2.7, strip.text = element_text(size=7)

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
  if (grepl(species_name, row_with_info[1]) == TRUE) {
    contig_name <- sub(paste0('_', species_name), '', row_with_info[1])
    contig_name <- sub('.p[1-9]*', '', contig_name)
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
  }
  return(new_df)
}

combined_data_transc <- NULL
for (row in 1:nrow(target_info)) {
  if (target_info[row, 'species'] == 'Ecy'){
    combined_data_transc <- 
      rbind(combined_data_transc,
            connect_counts_and_info(target_info[row,], counts_ecy, 'Ecy'))
  }  
  else if (target_info[row, 'species'] == 'Eve') {
    combined_data_transc <- 
      rbind(combined_data_transc,
            connect_counts_and_info(target_info[row,], counts_eve, 'Eve'))
  }
  else if (target_info[row, 'species'] == 'Gla') {
    combined_data_transc <- 
      rbind(combined_data_transc,
            connect_counts_and_info(target_info[row,], counts_gla, 'Gla'))
  }
}

colnames(combined_data_transc)[8] <- 'counts'
colnames(combined_data_transc)[9] <- 'sample'
combined_data_transc$condition <- ifelse(grepl('6C', combined_data_transc$sample), 
                                  '6C', '24C')

### 6. Plot the counts of the found target proteins:
to_plot <- subset(combined_data_transc, orthogroup == 'OG0000002')
to_plot$annotation <- sub(' isoform[^,]*,?', '', to_plot$annotation)
to_plot$annotation <- sub(' partial', '', to_plot$annotation)
to_plot$annotation <- sub('PREDICTED: ', '', to_plot$annotation)
to_plot$annotation <- sub('-like', '', to_plot$annotation)
to_plot$all_labels <- sprintf("%s %s|%s, %s", to_plot$species, to_plot$protein_group,
                              to_plot$annotation, to_plot$hsp70_type)

f <- function(x) {
  sapply(strsplit(x, '|', fixed=T), `[`, 2)
}
my_labeles <- c('Eve' = 'E. verrucosus', 
                'Ecy' = 'E. cyaneus', 
                'Gla' = 'G. lacustris')
my_colors <- c('#7570b3', '#1b9e77', '#d95f02')
ggplot(to_plot, aes(x = condition, y = counts)) +
  geom_boxplot(aes(color = species, fill = species), 
               outlier.alpha = 0, alpha = 0.4) +
  geom_jitter(aes(color = species)) +
  scale_color_manual('Species', values = my_colors,
                     labels = my_labeles) +
  scale_fill_manual('Species', values = my_colors,
                    labels = my_labeles) +
  facet_wrap(~all_labels, labeller = as_labeller(f), ncol = 3) +
  ylab('Scaled counts') +
  xlab('') +
  theme_bw() +
  theme(legend.text = element_text(face = 'italic', size = 11),
        axis.title.y = element_text(size = 12)
        ,
        strip.text = element_text(size=6))

dir_to_save <- '~/labeglo2/MS_results/390/withDBfromRNAspades/hsps_orthologes/hsps_transcriptomics/'
ggsave(file.path(dir_to_save, 'og7.png'),
       #scale = 1.2) 
       width = 9, height = 2.7)
ggsave(file.path(dir_to_save, 'og7.pdf'),
       #scale = 1.2)
       width = 9, height = 2.7)
# og0 - width = 10, height = 6
# og2 - width = 8, height = 4.5
# og1 - width = 8, height = 6
# og5 - width = 9, height = 2.7
# og7 - width = 9, height = 2.5, strip.text = element_text(size=6)
# og11, og13, og15, og16 -  width = 9, height = 2.7
# og19 - width = 6, height = 2.7, strip.text = element_text(size=7)

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
combined_data <- rbind(combined_data_proteomics, combined_data_transc)
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
}
combined_data <- unique(combined_data)

to_plot <- subset(combined_data, orthogroup == "OG0000002")
to_plot$annotation <- sub(' isoform[^,]*,?', '', to_plot$annotation)
to_plot$annotation <- sub(' partial', '', to_plot$annotation)
to_plot$annotation <- sub('PREDICTED: ', '', to_plot$annotation)
to_plot$annotation <- sub('-like', '', to_plot$annotation)
to_plot$all_labels <- sprintf("%s%s|%s, %s",
                              to_plot$species, to_plot$protein_group,
                              to_plot$annotation, to_plot$hsp70_type)

to_plot_ecy <- subset(to_plot, species == 'Ecy')
to_plot_eve <- subset(to_plot, species == 'Eve')
to_plot_gla <- subset(to_plot, species == 'Gla')

plot_target_proteins <- function(to_plot){
  f <- function(x) {
    sapply(strsplit(x, '|', fixed=T), `[`, 1)
  }
  ggplot(to_plot, aes(x = method, y = values, color = condition)) +
  geom_point(position=position_jitterdodge(dodge.width=.8)) +
  geom_boxplot(aes(fill = condition), outlier.alpha = 0, alpha = 0.4) +
  facet_wrap(~all_labels, labeller = as_labeller(f), nrow = 1) +
  theme_bw()
}

(ecy_plot <- plot_target_proteins(to_plot_ecy))
(eve_plot <- plot_target_proteins(to_plot_eve))
(gla_plot <- plot_target_proteins(to_plot_gla))

ggarrange(eve_plot, gla_plot, ecy_plot, nrow = 3, widths = c(5, 4 ,1)) 

library(cowplot)
plot_grid(eve_plot, gla_plot, ecy_plot, ncol = 1, rel_widths = c(5, 4, 1))

ggdraw() +
  draw_plot(eve_plot, 0, .66, 1, .33) +
  draw_plot(gla_plot, 0, .33, 0.85, .33) +
  draw_plot(ecy_plot, 0, 0, 0.3, .33)

library(ggh4x)
f <- function(x) {
  sapply(strsplit(x, '|', fixed=T), `[`, 2)
}
to_plot$species <- factor(to_plot$species, levels = c('Eve', 'Gla', 'Ecy'), 
                          labels = c('E.verrucosus', 'G.lacustris', 
                                     'E.cyaneus'))
to_plot$species_italic <- sprintf('italic(%s)', to_plot$species)
to_plot$species_italic <- factor(to_plot$species_italic, 
                                 levels = c('italic(E.verrucosus)', 
                                            'italic(G.lacustris)', 
                                            'italic(E.cyaneus)'))
  
ggplot(to_plot, aes(x = method, y = values, color = condition)) +
  geom_point(position=position_jitterdodge(dodge.width=.8)) +
  geom_boxplot(aes(fill = condition), outlier.alpha = 0, alpha = 0.4) +
  facet_nested_wrap(species_italic ~ all_labels, 
                    labeller = labeller(all_labels = as_labeller(f), 
                                        species_italic = label_parsed), 
                    nrow = 2, scales = 'free') +
  scale_color_manual('Condition:', values = c('#ca0020', '#0571b0')) +
  scale_fill_manual('Condition:', values = c('#ca0020', '#0571b0')) +
  xlab('') +
  ylab('Scaled absolute values') +
  theme_bw()

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

combined_data$pg_unique_name <- paste0(combined_data$protein_group, '_', 
                                       combined_data$species)
combined_data$contig_unique_name <- sub('\\.p[1-9]*_...$', '', 
                                        combined_data$protein)
combined_data$contig_unique_name <- paste0(combined_data$contig_unique_name, '_', 
                                       combined_data$species)

combined_data <- cbind(combined_data,
                       total_proteins[match(combined_data$pg_unique_name, total_proteins$unique_name),][6])
combined_data <- cbind(combined_data, 
                       total_transc[match(combined_data$contig_unique_name, total_transc$unique_name),][6])

total_proteins[total_proteins$protein == 'NODE_9313_length_2618_cov_250.297392_g5535_i0.p1_Gla',]

combined_data$sign2plot <- 
  ifelse(combined_data$FDR_recalc < 0.05 & combined_data$padj < 0.05,
         '< 0.05', 
         ifelse(combined_data$padj > 0.05 & combined_data$method == 'RNAseq',
                '> 0.05', 
                ifelse(combined_data$FDR_recalc > 0.05 & combined_data$method == 'MS/MS',
                                 '> 0.05', '< 0.05')))

to_plot <- subset(combined_data, orthogroup == "OG0000002")
to_plot$annotation <- sub(' isoform[^,]*,?', '', to_plot$annotation)
to_plot$annotation <- sub(' partial', '', to_plot$annotation)
to_plot$annotation <- sub('PREDICTED: ', '', to_plot$annotation)
to_plot$annotation <- sub('-like', '', to_plot$annotation)
to_plot$all_labels <- sprintf("%s%s|%s, %s",
                              to_plot$species, to_plot$protein_group,
                              to_plot$annotation, to_plot$hsp70_type)
f <- function(x) {
  sapply(strsplit(x, '|', fixed=T), `[`, 2)
}
to_plot$species <- factor(to_plot$species, levels = c('Eve', 'Gla', 'Ecy'), 
                          labels = c('E.verrucosus', 'G.lacustris', 
                                     'E.cyaneus'))
to_plot$species_italic <- sprintf('italic(%s)', to_plot$species)
to_plot$species_italic <- factor(to_plot$species_italic, 
                                 levels = c('italic(E.verrucosus)', 
                                            'italic(G.lacustris)', 
                                            'italic(E.cyaneus)'))

to_plot <- na.omit(to_plot)
to_plot$condition <- factor(to_plot$condition, levels = c('24C', '6C'),
                            labels = c('24 °C', '6 °C'))
to_plot$species_italic <- factor(to_plot$species_italic,
                                 levels = c('italic(E.cyaneus)',
                                            'italic(E.verrucosus)',
                                            'italic(G.lacustris)'))
##### if do not want ', *' in the end of annotations:
to_plot$all_labels <- sub(', \\*', '', to_plot$all_labels)
#####
ggplot(to_plot, aes(x = method, y = values, color = condition)) +
  geom_point(position=position_jitterdodge(dodge.width=.8)) +
  geom_boxplot(aes(fill = condition,
                   linetype = sign2plot), outlier.alpha = 0, alpha = 0.4) +
  facet_nested_wrap(species_italic ~ all_labels, 
                    labeller = labeller(all_labels = as_labeller(f), 
                                        species_italic = label_parsed), 
                    scales = 'free', ncol = 4) +
  scale_color_manual('Condition:', values = c('#ca0020', '#0571b0')) +
  scale_fill_manual('Condition:', values = c('#ca0020', '#0571b0')) +
  scale_linetype('adj. p-value:') +
  xlab('') +
  ylab('Scaled absolute values') +
  theme_bw() 
#  theme(strip.text = element_text(size=5.5))

dir_to_save <- '/home/polina/labeglo2/MS_results/390/withDBfromRNAspades/hsps_orthologes'
ggsave(file.path(dir_to_save, 'og2_proteinsWITHtranscripts.png'),
       #scale = 1.2) 
       width = 10, height = 5)
ggsave(file.path(dir_to_save, 'og2_proteinsWITHtranscripts.pdf'),
       #scale = 1.2)
       width = 10.2, height = 5)
# og0, og1 - width = 10, height = 5
# og2 - width = 8, height = 5
# og11, og5, og13, og16, og15 - width = 8, height = 3

