library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(rstatix)

dir <- '~/labeglo2/proteome_transcr_comparision'
transcr <- 
  read.csv('~/labeglo2/proteome_transcr_comparision/Ecy_transcr_24vs6_all.csv', 
           sep = '\t')
proteins <- 
  read.csv('~/labeglo2/proteome_transcr_comparision/Ecy_AllProteins_24vs6_proteinGroups_separatelyAnalyzed.csv',
          sep = '\t', header = T)

annot <- read.csv(file.path(dir, 'contigs_whole_annot_Ecy.csv'), sep = '\t')
#colnames(proteins) <- c('protein', 'geneSymbol', 'logFC', 'pvalue', 'FDR')

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

joined <- inner_join(transcr, proteins_long, by = c('contig' = 'protein_clip'))

joined_clip <- joined[c(1, 3, 5, 7, 8, 9, 11, 12, 13)]

filtered_joined_clip <- filter(joined_clip, abs(log2FoldChange) > 2 | abs(logFC) > 0.5)
#filtered_joined_clip

set.seed(333)
ggplot(joined_clip, aes(log2FoldChange, logFC)) +
  geom_point(alpha=.5, color='gray70') +
  geom_point(data=filtered_joined_clip, color='red') +
  geom_text_repel(aes(label= geneSymbol),
            data = filtered_joined_clip) +
  xlab('log2FC (24°C/6°C) transcriptome') +
  ylab('log2FC (24°C/6°C) proteome') +
  theme_light() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))

#ggsave(file.path(dir, 'transcr_proteome_logfc.png'), scale = 3)

# it seems that I have a lot of duplicated proteins 
# (no wonder - protein groups were splitted and transcripts matched to 
# different proteins from the group which is basically the same protein)

absmax <- function(x) { x[which.max(abs(x))][1] }
joined_clip_merged <- joined_clip %>% 
  group_by(protein_group) %>%
  mutate(best_tlfc = absmax(log2FoldChange)) %>% # tlfc - transcriptome lfc
  mutate(pvalue_tran = min(pvalue.x)) %>%
  mutate(transcr_groups = paste0(contig, collapse = ';')) %>%
  mutate(proteome_groups = paste0(protein, collapse = ';')) %>%
  mutate(sign = ifelse(
    pvalue.y < 0.05 & !is.na(pvalue.y) & pvalue_tran < 0.05 & !is.na(pvalue_tran), 
                       "< 0.05 (both)",
      ifelse(pvalue_tran < 0.05 & !is.na(pvalue_tran), "< 0.05 (transcriptome)", 
             ifelse(pvalue.y < 0.05 & !is.na(pvalue.y), 
                    "< 0.05 (proteome)", "> 0.05 (both)")))) %>%
  ungroup() %>%
  dplyr::select(!c(contig, log2FoldChange, protein, pvalue.x)) %>%
  unique()

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

cor_test_res <- cor.test(joined_clip_merged$best_tlfc, joined_clip_merged$logFC)

ggplot(joined_clip_merged, aes(best_tlfc, logFC)) +
  geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
  geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
  geom_point(aes(shape = sign), alpha=.7, color='gray70') +
  geom_smooth(method = 'lm') +
  annotate(geom='text', x = -10, y = 1.8, hjust = 0, size = 4.5,
           label = paste0('r2 = ', round(cor_test_res$estimate, 4), '\n',
                          'p-value ', p_format(cor_test_res$p.value,
                                               accuracy = 0.0001))) +
#  geom_point(aes(shape = sign), data=filtered_joined_clip_merged, color='red') +
  scale_shape_manual('p-value:', values = c(8, 16, 17, 1)) +
  xlab('log2FC (24°C/6°C) transcriptome') +
  ylab('log2FC (24°C/6°C) proteome') +
  theme_light() +
  guides(shape = guide_legend(override.aes = list(size=4))) +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15))

ggsave(file.path(dir, 'transcr_proteome_logfc_Ecy_pvalues_24Cvs6C_cor.png'),
       scale = 1.5)

# visualize all transcripts versus all proteins 
# (so, you can see matched proteins/transcripts and not-matched ones)

joined_full <- full_join(transcr, proteins_long, by = c('contig' = 'protein_clip'))

joined_full$tag <- ifelse(is.na(joined_full$protein), 'transcriptome_only',
                          ifelse(is.na(joined_full$baseMean), 'proteome_only',
                          'matched'))

joined_full$logFC_transc <- replace_na(joined_full$log2FoldChange, 0)
joined_full$logFC_proteome <- replace_na(joined_full$logFC, 0)

filtered_joined_full <- 
  filter(joined_full,
         (logFC_transc >5 | logFC_transc < -5 & tag == 'transcriptome_only')  |
           (abs(logFC_proteome) > 0.8 & tag == 'proteome_only') )

ixs <- is.na(filtered_joined_full$protein_group)
filtered_joined_full[ixs, 'protein_group'] <- 
  seq(max(filtered_joined_full$protein_group, na.rm = T) * 100, 
      max(filtered_joined_full$protein_group, na.rm = T) * 100 + sum(ixs) - 1)

filtered_joined_full <- filtered_joined_full %>%
  group_by(protein_group) %>% slice_head(n = 1)

ixs <- is.na(filtered_joined_full$geneSymbol)
filtered_joined_full[ixs, 'geneSymbol'] <- annot[match(filtered_joined_full$contig[ixs],
                                                       annot$contig),]$annot
filtered_joined_full$geneSymbol <- sub('PREDICTED: |', '', 
                                       filtered_joined_full$geneSymbol)
filtered_joined_full$geneSymbol <- sub('-like.*', '', 
                                       filtered_joined_full$geneSymbol)
set.seed(356)
ggplot(joined_full, aes(logFC_transc, logFC_proteome)) +
  geom_point(aes(color = tag), alpha=.3) +
  scale_color_manual('Occurance', values = c('#354E6C', '#E64241', '#EAAC31')) +
  theme_light() +
  geom_text_repel(aes(label = ifelse(grepl('uncharacterized|hypothetical|\\*',
                                           geneSymbol),
                                     '', geneSymbol)),
                  data = filtered_joined_full, 
                  segment.colour = 'grey50', 
                  max.overlaps = Inf) +
  xlab('log2FC (24°C/6°C) transcriptome') +
  ylab('log2FC (24°C/6°C) proteome') +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15))

ggsave(file.path(dir, 'transcr_proteome_logfc_Gla_allObservations_24vs6Cafter.png'), 
       scale = 3)
