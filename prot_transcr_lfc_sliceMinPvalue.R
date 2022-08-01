### Heat shock experiment analysis ###
#### Merge the proteomic and transcriptomic DA/DE analysis data and plot figures ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(rstatix)

species <- 'Eve'
dir <- '~/labeglo2/proteome_transcr_comparision'
dir <- '~/labeglo2/proteome_transcr_comparision/3h/'
transcr <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/', species,
                  '_transcr_24vs6_all.csv'), 
           sep = '\t') # 24h
transcr <-
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/3h/', species,
                  '_transcr_24vs6_all_3h.csv'), 
           sep = '\t') # 3h
transcr <- subset(transcr, pvalue < 0.05)
proteins <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/', species,
                  '_AllProteins_24vs6_proteinGroups_separatelyAnalyzed.csv'),
           sep = '\t', header = T) # 24h
proteins <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/3h/', species,
                  '_AllProteins_24vs6_proteinGroups_separatelyAnalyzed_3h.csv'),
           sep = '\t', header = T) # 3h
proteins <- subset(proteins, pvalue < 0.05)

annot <- read.csv(file.path(dir, 
                            paste0('contigs_whole_annot_', 
                                   species,'.csv')), sep = '\t') # 24h

annot <- read.csv(file.path(dir, 
                            paste0('contigs_whole_annot_', 
                                   species,'_3h.csv')), sep = '\t') # 3h

proteins$protein_group <- 1:nrow(proteins)
proteins_sep <- separate(proteins, protein, ';', into=paste0('key', 1:31), 
                         fill='right')
proteins_long <- pivot_longer(proteins_sep, starts_with('key'),
                              names_to = NULL, values_to = 'protein') %>%
  filter(!is.na(protein))

proteins_long$protein_clip <- sub('\\.p[0-9]+_...$', '', proteins_long$protein)
length(unique(proteins_long$protein_clip))
proteins_long <- group_by(proteins_long, protein_clip) %>%
  slice_min(pvalue, n = 1, with_ties = F)

joined <- inner_join(transcr, proteins_long, by = c('contig' = 'protein_clip'))

joined_clip <- joined[c(1, 3, 5, 7, 8, 9, 11, 12, 13)]

filtered_joined_clip <- filter(joined_clip, abs(log2FoldChange) > 2 | abs(logFC) > 0.5)

##### add ORF length to all complete ORFs
#orf_type <- read.table(file = file.path(dir, 'all_species_orf_type.tsv'))
#colnames(orf_type) <- c('protein_name', 'orf', 'length')
#complete_orf <- subset(orf_type, orf == 'complete')
#complete_orf[match(proteins_long$protein, complete_orf$protein_name),]
#proteins_long_complete <- subset(proteins_long, protein %in% complete_orf$protein_name)
#proteins_long_complete$orf_length <- 
#complete_orf[match(proteins_long_complete$protein, complete_orf$protein_name),]$length

#joined <- inner_join(transcr, proteins_long_complete, by = c('contig' = 'protein_clip'))

#joined_clip <- joined[c(1, 3, 5, 7, 8, 9, 11, 12, 13, 14)]

# it seems that I have a lot of duplicated proteins 
# (no wonder - protein groups were splitted and transcripts matched to 
# different proteins from the group which is basically the same protein)

absmax <- function(x) { x[which.max(abs(x))][1] }
joined_clip_merged <- joined_clip %>% 
  group_by(protein_group) %>%
  mutate(transcr_groups = paste0(contig, collapse = ';')) %>%
  mutate(proteome_groups = paste0(protein, collapse = ';')) %>%
  slice_min(pvalue.x, n = 1, with_ties = F) %>%
  mutate(best_tlfc = log2FoldChange) %>% # tlfc - transcriptome lfc
  mutate(pvalue_tran = pvalue.x) %>%
  mutate(sign = ifelse(
    pvalue.y < 0.05 & !is.na(pvalue.y) & pvalue_tran < 0.05 & !is.na(pvalue_tran), 
    "< 0.05 (both)",
    ifelse(pvalue_tran < 0.05 & !is.na(pvalue_tran), "< 0.05 (transcriptome)", 
           ifelse(pvalue.y < 0.05 & !is.na(pvalue.y), 
                  "< 0.05 (proteome)", "> 0.05 (both)")))) %>%
  ungroup() %>%
  dplyr::select(!c(contig, log2FoldChange, protein, pvalue.x))

joined_clip_merged$sign <- factor(joined_clip_merged$sign, 
                                  levels = c("< 0.05 (both)", 
                                             "< 0.05 (proteome)", 
                                             "< 0.05 (transcriptome)", 
                                             "> 0.05 (both)"))

filtered_joined_clip_merged <- 
  filter(joined_clip_merged, 
         abs(best_tlfc) > 2.5 | abs(logFC) > 0.4 ) # Ecy

filtered_joined_clip_merged <- 
  filter(joined_clip_merged, 
         abs(best_tlfc) > 2 | abs(logFC) > 0.5 ) # Gla

filtered_joined_clip_merged <- 
  filter(joined_clip_merged, 
         abs(best_tlfc) > 5 | abs(logFC) > 0.5 ) # Eve

filtered_joined_clip_merged$geneSymbol <- sub('PREDICTED: ', '', 
                                              filtered_joined_clip_merged$geneSymbol)
filtered_joined_clip_merged$geneSymbol <- sub('-like.*', '', 
                                              filtered_joined_clip_merged$geneSymbol)

set.seed(365)
ggplot(joined_clip_merged, aes(best_tlfc, logFC)) +
  geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
  geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
  geom_point(aes(shape = sign), alpha=.5, color='gray70') +
  geom_text_repel(aes(label = ifelse(grepl('uncharacterized|hypothetical|\\*',
                                           geneSymbol),
                                     '', geneSymbol)),
                  data = filtered_joined_clip_merged, 
                  segment.colour = 'grey50',
                  max.overlaps = Inf, size = 4) +
  geom_point(aes(shape = sign), data=filtered_joined_clip_merged, color='red') +
  scale_shape_manual('p-value:', values = c(8, 16, 17, 1)) +
  xlab('log2FC (24°C/6°C) transcriptome') +
  ylab('log2FC (24°C/6°C) proteome') +
  theme_light() +
  guides(shape = guide_legend(override.aes = list(size=4))) +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15))

ggsave(file.path(dir, 'transcr_proteome_logfc_Gla_pvalues_24Cvs6Cafter.png'), 
       scale = 3)

# figure with correlation curve and cor_test results

joined_clip_merged_ <- joined_clip_merged
#joined_clip_merged <- joined_clip_merged_
#joined_clip_merged <- subset(joined_clip_merged_, orf_length >= 200)
cor_test_res <- cor.test(joined_clip_merged$best_tlfc, joined_clip_merged$logFC)

ggplot(joined_clip_merged, aes(best_tlfc, logFC)) +
  geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
  geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
  geom_point(aes(shape = sign), alpha=.7, color='gray70') +
  geom_smooth(method = 'lm') +
  annotate(geom='text', x = -5.5, y = 1.6, hjust = 0, size = 6,
           label = paste0('r2 = ', round(cor_test_res$estimate, 4), '\n',
                          'p-value ', p_format(cor_test_res$p.value,
                                               accuracy = 0.00001))) +
  #  geom_point(aes(shape = sign), data=filtered_joined_clip_merged, color='red') +
  scale_shape_manual('p-value:', values = c(8, 16, 17, 1)) +
  #xlab('log2FC (24°C/6°C) transcriptome') + # 24h
  xlab('log2FC (24°C/6°C) transcriptome, 3 hours exposure') + # 3h
  #ylab('log2FC (24°C/6°C) proteome') +
  ylab('log2FC (24°C/6°C) proteome, 24 hours exposure') +
  theme_light() +
  guides(shape = guide_legend(override.aes = list(size=4))) +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15))

ggsave(file.path(dir, paste0('transcr_proteome_logfc_', 
                             species,'_pvalues_24Cvs6C_cor_3h.png')),
       scale = 1.9)

#write.table(joined_clip_merged, 
#            file = file.path(dir, paste0(species, '_3h_table_for_cor_plot_All_sliceMinPValue.csv')),
#            sep = '\t')

### To draw plots only those proteins and transcripts that have p-value < 0.05
cor_test_res <- cor.test(joined_clip_merged$best_tlfc, joined_clip_merged$logFC)

joined_clip_merged$geneSymbol <- sub('PREDICTED: |', '', 
                                     joined_clip_merged$geneSymbol)
joined_clip_merged$geneSymbol <- sub('-like.*', '', 
                                     joined_clip_merged$geneSymbol)
joined_clip_merged$geneSymbol <- sub('LOW QUALITY PROTEIN: ', '', 
                                     joined_clip_merged$geneSymbol)

var_width <- 25
joined_clip_merged <- mutate(joined_clip_merged, 
                             pretty_varname = stringr::str_wrap(joined_clip_merged$geneSymbol, 
                                                                width = var_width))

ggplot(joined_clip_merged, aes(best_tlfc, logFC)) +
  geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
  geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
  geom_point(color='#0099E5') +
  geom_smooth(method = 'lm', color = 'grey65', fill = 'grey85') +
  annotate(geom='text', x = -1.5, y = 1.6, hjust = 0, size = 7,
           label = paste0('r2 = ', round(cor_test_res$estimate, 4), '\n',
                          'p-value ', p_format(cor_test_res$p.value,
                                               accuracy = 0.01))) +
  geom_text_repel(aes(segment.color = 'grey70',
                      label = ifelse(grepl('uncharacterized|hypothetical|unknown|\\*',
                                           pretty_varname),
                                     '', pretty_varname)),
                  max.overlaps = Inf, size = 4, max.time = 1, box.padding = 0.7,
                  seed = 546) +
  #xlab('log2FC (24°C/6°C) transcriptome') + # 24h
  xlab('log2FC (24°C/6°C) transcriptome, 3 hours exposure') + # 3h
  #ylab('log2FC (24°C/6°C) proteome') +
  ylab('log2FC (24°C/6°C) proteome, 24 hours exposure') +
  theme_light() +
  theme(axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text = element_text(size = 14))

write.table(joined_clip_merged, 
            file = file.path(dir, paste0(species, 
                             '_3h_table_for_cor_plot_onlySign_sliceMinPValue.csv')),
            sep = '\t')
