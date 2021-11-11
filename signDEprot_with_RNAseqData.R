library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(rstatix)
library(stringr)

species <- 'Eve'
dir <- '~/labeglo2/proteome_transcr_comparision'
dir <- '~/labeglo2/proteome_transcr_comparision/3h/'
transcr <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/', species, '_transcr_24vs6_all.csv'), 
           sep = '\t')
transcr <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/3h/', species, '_transcr_24vs6_all_3h.csv'), 
           sep = '\t') # 3h
#transcr <- subset(transcr, padj < 0.05)
proteins <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/', species, '_AllProteins_24vs6after_proteinGroups_separatelyAnalyzed.csv'),
           sep = '\t', header = T)
proteins <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/3h/', species, '_AllProteins_24vs6_proteinGroups_separatelyAnalyzed_3h.csv'),
           sep = '\t', header = T) # 3h
proteins <- subset(proteins, FDR_recalc < 0.05)

annot <- read.csv(file.path(dir, paste0('contigs_whole_annot_', species,'.csv')), 
                  sep = '\t')
annot <- read.csv(file.path(dir, paste0('contigs_whole_annot_', species,'_3h.csv')), 
                  sep = '\t') # 3h

proteins$protein_group <- 1:nrow(proteins)
proteins_sep <- separate(proteins, protein, ';', into=paste0('key', 1:30), 
                         fill='right')
proteins_long <- pivot_longer(proteins_sep, starts_with('key'),
                              names_to = NULL, values_to = 'protein') %>%
  filter(!is.na(protein))

proteins_long$protein_clip <- sub('\\.p[0-9]+_...$', '', proteins_long$protein)
length(unique(proteins_long$protein_clip))
proteins_long <- group_by(proteins_long, protein_clip) %>%
  slice_min(pvalue, n = 1, with_ties = F)

# here enrich_DE_forProteomeTranscriptomeComp.R can be run
joined <- inner_join(transcr, proteins_long, by = c('contig' = 'protein_clip'))

joined_clip <- joined[c(1, 3, 6, 7, 8, 9, 11, 12, 13)]

# it seems that I have a lot of duplicated proteins 
# (no wonder - protein groups were splitted and transcripts matched to 
# different proteins from the group which is basically the same protein)

absmax <- function(x) { x[which.max(abs(x))][1] }
joined_clip_merged <- joined_clip %>% 
  group_by(protein_group) %>%
#  mutate(best_tlfc = absmax(log2FoldChange)) %>% # tlfc - transcriptome lfc
  slice_max(abs(log2FoldChange), n = 1, with_ties = F) %>%
  mutate(padj_tran = min(padj)) %>%
  mutate(transcr_groups = paste0(contig, collapse = ';')) %>%
  mutate(proteome_groups = paste0(protein, collapse = ';')) %>%
  mutate(sign = ifelse(
    FDR_recalc < 0.05 & !is.na(FDR_recalc) & padj_tran < 0.05 & !is.na(FDR_recalc), 
    "< 0.05 (both)",
    ifelse(padj_tran < 0.05 & !is.na(padj_tran), "< 0.05 (transcriptome)", 
           ifelse(FDR_recalc < 0.05 & !is.na(FDR_recalc), 
                  "< 0.05 (proteome)", "> 0.05 (both)")))) %>%
  ungroup() %>%
#  dplyr::select(!c(contig, log2FoldChange, protein, padj)) %>%
  dplyr::select(!c(protein, padj)) %>%
  unique()

### Connect joined data with intensities (proteome) and counts (transcriptome) data
# 1. Intensities
intensities <- 
  read.table(paste0('~/labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_', 
                    tolower(species),'/intensities_after_slNorm_', 
                    tolower(species),'.csv'), 
            header = T)
intensities <- 
  read.table(paste0('~/labeglo2/MS_results/390/3h/', species,'/intensities_after_slNorm_', 
                    tolower(species),'.csv'), 
             header = T) # 3h
data_prep <- function(data_intensities){
  row.names(data_intensities) <- data_intensities$protein_group
  data_intensities <- data_intensities[-1]
  data_intensities_scaled <- apply(data_intensities, 1, scale)
  data_intensities_scaled <- t(data_intensities_scaled)
  colnames(data_intensities_scaled) <- colnames(data_intensities)
  return(data_intensities_scaled)
}

intensities <- data_prep(intensities)

# 2. Counts:
count_dir <- paste0('~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/', 
              species, '/counts') # specify path to your samples
count_dir <- paste0('~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/3h/', 
                    species, '/count')

countFiles <- list.files(count_dir, full.names = T)
countFiles
counts <- lapply(countFiles, function(countsFile) {
  f <- file.path(countsFile, 'quant.sf')
  read.table(f, sep="\t", header=1, row.names = 1, 
             stringsAsFactors = F, comment.char = "#")
})
counts <- lapply(counts, function(countsTable) countsTable[, 4, drop=F])

counts <- do.call(cbind, counts)

if (species == 'Gla'){
  sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3',
                    '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4')
} else {
  sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3', '24C_rep4',
                    '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4')
}
if (species == 'Eve'){
  sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3','24C_rep4',
                    '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4')
} else {
  sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3',
                    '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4')
} # 3h
colnames(counts) <- sample_names

data_prep_counts <- function(data_counts){
  data_counts_scaled <- apply(data_counts, 1, scale)
  data_counts_scaled <- t(data_counts_scaled)
  colnames(data_counts_scaled) <- colnames(data_counts)
  return(data_counts_scaled)
}
counts <- data_prep_counts(counts)

# 3. Connect the data:

connect_info <- function(row_with_info, intens_data, counts_data){
  add_new_rows <- function(list_with_values, df2update, info_row){
    for (el in 1:length(list_with_values)){
      if (length(list_with_values) == 0) {
        break
      }
      new_row <- info_row
      new_row[[length(new_row) + 1]] <- list_with_values[el]
      new_row[[length(new_row) + 1]] <- names(list_with_values[el])
      df2update <- rbind(df2update, new_row)
    }
    return(df2update)
  }
  new_df <- NULL
  pg_name <- row_with_info[12]
  contig_name <- row_with_info[1]
  intens_values <- intens_data[row.names(intens_data) == as.character(pg_name),]
  count_values <- counts_data[row.names(counts_data) == as.character(contig_name),]
  new_df <- add_new_rows(intens_values, new_df, row_with_info)
  new_df <- add_new_rows(count_values, new_df, row_with_info)
  return(new_df)
}
# joined_clip_merged df does not contain full protein_group name, to retrieve:
retrieved_pg <- match(joined_clip_merged$protein_group, proteins$protein_group)
joined_clip_merged$pg_name <- proteins$protein[retrieved_pg]
#
combined_data <- NULL
for (row in 1:nrow(joined_clip_merged)){
  combined_data <- rbind(combined_data, connect_info(joined_clip_merged[row,],
                                                     intensities,
                                                     counts))
}
colnames(combined_data)[13] <- 'values'
colnames(combined_data)[14] <- 'sample'
combined_data$condition <- 
  ifelse(grepl('6C|CK1|CK2|LK1|LK2|VrK1|VrK2', combined_data$sample), 
         '6C', '24C')

combined_data$method <- ifelse(grepl('rep', combined_data$sample),
                               'RNAseq', 'MS/MS')
combined_data$geneSymbol <- sub('PREDICTED: ', '', combined_data$geneSymbol)
combined_data$geneSymbol <- sub('-like', '', combined_data$geneSymbol)
combined_data$geneSymbol <- sub(', partial', '', combined_data$geneSymbol)
var_width <- 25
#combined_data <- tmp
combined_data <- mutate(combined_data, 
                        pretty_varname = str_wrap(combined_data$geneSymbol, 
                                                  width = var_width))
combined_data$all_labels <- sprintf('%s|%s', combined_data$pg_name, 
                                    combined_data$pretty_varname)
combined_data$all_labels2 <- factor(combined_data$all_labels,
                             levels=unique(
                             combined_data$all_labels[order(combined_data$pretty_varname)]))

combined_data$sign2plot <- 
  ifelse(combined_data$sign == "< 0.05 (both)",
        '< 0.05', 
         ifelse(combined_data$sign == "< 0.05 (proteome)" & combined_data$method == 'RNAseq',
                '> 0.05', '< 0.05'))

combined_data$method <- factor(combined_data$method, 
                               levels = c('RNAseq', 'MS/MS'))
combined_data$condition <- factor(combined_data$condition, 
                                  levels = c('6C', '24C'),
                                  labels = c('6 °C', '24.6 °C'))
combined_data_up <- subset(combined_data, logFC > 0)
write.table(combined_data_up, 
            file = file.path(dir, 
                   paste(species, '_AllupProteins_joinedWithTranscripts.csv')))
combined_data_down <- subset(combined_data, logFC < 0)
write.table(combined_data_down, 
            file = file.path(dir, 
            paste(species, '_AlldownProteins_joinedWithTranscripts.csv')))

f <- function(x) {
  sapply(strsplit(x, '|', fixed=T), `[`, 2)
}
ggplot(combined_data, aes(x = method, y = values, color = condition)) +
  geom_point(position=position_jitterdodge(dodge.width=1),
             size = 0.7) +
  geom_boxplot(aes(fill = condition, 
                   linetype = sign2plot), 
               outlier.alpha = 0, alpha = 0.4) +
  facet_wrap(~all_labels2, labeller = as_labeller(f), ncol = 5) +
  scale_linetype('adj. p-value:') +
  scale_color_manual('Condition:', values = c('#0571b0', '#ca0020')) +
  scale_fill_manual('Condition:', values = c('#0571b0', '#ca0020')) +
  theme_bw() +
  ylab('Scaled absolute values') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 8, margin = margin(2,2,2,2)))

dir_to_save <- '~/labeglo2/proteome_transcr_comparision/'
ggsave(file.path(dir_to_save, paste0(tolower(species), 
                                     '_DEproteinsWITHtranscripts_up.png')),
       scale = 0.8,
       width = 11.5, height = 7)
ggsave(file.path(dir_to_save, paste0(tolower(species), 
                                     '_DEproteinsWITHtranscripts_up.pdf')),
       scale = 0.8,
       width = 12, height = 7)
# Ecy: width = 12.5, height = 7
# Gla: width = 9.5, height = 5
# Eve: width = 12.5, height = 7