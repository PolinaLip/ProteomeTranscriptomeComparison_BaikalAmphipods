### Heat shock experiment analysis ###
#### Merge the proteomic and transcriptomic DA/DE analysis data and plot figures ####
## Here is from the perspective of DE transcripts ##
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
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/', species,
                  '_transcr_24vs6_all.csv'), 
           sep = '\t') # 24h
transcr <- 
  read.csv(file.path(dir,'Eve_wo_EveB24_4', 'Eve_transcr_24vs6_all_24h_wo__EveB24_4.csv'), 
           sep = '\t')
transcr <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/3h/', species, '_transcr_24vs6_all_3h.csv'), 
           sep = '\t') # 3h
#transcr <- subset(transcr, pvalue < 0.05)
proteins <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/', species,
                  '_AllProteins_24vs6_proteinGroups_separatelyAnalyzed.csv'),
           sep = '\t', header = T) # 24h
proteins <- 
  read.csv(paste0('~/labeglo2/proteome_transcr_comparision/3h/', species, '_AllProteins_24vs6_proteinGroups_separatelyAnalyzed_3h.csv'),
           sep = '\t', header = T) # 3h
#proteins <- subset(proteins, pvalue < 0.05)

annot <- read.csv(file.path(dir, 
                            paste0('contigs_whole_annot_', 
                                   species,'.csv')), sep = '\t') # 24h
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

joined <- inner_join(transcr, proteins_long, by = c('contig' = 'protein_clip'))

joined_clip <- joined[c(1, 3, 6, 7, 8, 9, 11, 12, 13)]

filtered_joined_clip <- filter(joined_clip, abs(log2FoldChange) > 2 | abs(logFC) > 0.5)

# it seems that I have a lot of duplicated proteins 
# (no wonder - protein groups were splitted and transcripts matched to 
# different proteins from the group which is basically the same protein)

absmax <- function(x) { x[which.max(abs(x))][1] }
joined_clip_merged <- (joined_clip %>% 
  group_by(protein_group) %>%
  #mutate(best_tlfc = absmax(log2FoldChange)) %>% # tlfc - transcriptome lfc
  mutate(pvalue_tran = min(padj)) %>%
  mutate(transcr_groups = paste0(contig, collapse = ';')) %>%
  mutate(proteome_groups = paste0(protein, collapse = ';')) %>%
  slice_min(padj, n = 1, with_ties = F) %>%
  #slice_max(abs(log2FoldChange), n = 1, with_ties = F) %>%
  mutate(sign = ifelse(
    FDR_recalc < 0.05 & !is.na(FDR_recalc) & pvalue_tran < 0.05 & !is.na(pvalue_tran), 
    "< 0.05 (both)",
    ifelse(pvalue_tran < 0.05 & !is.na(pvalue_tran), "< 0.05 (transcriptome)", 
           ifelse(FDR_recalc < 0.05 & !is.na(FDR_recalc), 
                  "< 0.05 (proteome)", "> 0.05 (both)")))) %>%
  ungroup() %>%
  dplyr::select(!c(protein, padj)))

# pvalue.y is a nonadjusted pvalue from proteomic dataset

joined_clip_merged$sign <- factor(joined_clip_merged$sign, 
                                  levels = c("< 0.05 (both)", 
                                             "< 0.05 (proteome)", 
                                             "< 0.05 (transcriptome)", 
                                             "> 0.05 (both)"))

filtered_joined_clip_merged <- 
  filter(joined_clip_merged, 
         abs(log2FoldChange) > 2 | abs(logFC) > 0.4 ) # Ecy, logFC - logFC from proteomic dataset

filtered_joined_clip_merged <- 
  filter(joined_clip_merged, 
         abs(log2FoldChange) > 2 | abs(logFC) > 0.5 ) # Gla

filtered_joined_clip_merged <- 
  filter(joined_clip_merged, 
         abs(log2FoldChange) > 5 | abs(logFC) > 0.5 ) # Eve

filtered_joined_clip_merged$geneSymbol <- sub('PREDICTED: ', '', 
                                              filtered_joined_clip_merged$geneSymbol)
filtered_joined_clip_merged$geneSymbol <- sub('-like.*', '', 
                                              filtered_joined_clip_merged$geneSymbol)

ggplot(joined_clip_merged, aes(log2FoldChange, logFC)) +
  geom_hline(yintercept = 0, color = '#387490', alpha = .7) +
  geom_vline(xintercept = 0, color = '#387490', alpha = .7) +
  geom_point(aes(shape = sign), alpha=.5, color='gray70') +
  geom_text_repel(aes(label = ifelse(grepl('uncharacterized|hypothetical|\\*',
                                           geneSymbol),
                                     '', geneSymbol)),
                  data = filtered_joined_clip_merged, 
                  segment.colour = 'grey50',
                  max.overlaps = Inf, size = 3, min.segment.length = 0, seed = 88) +
  geom_point(aes(shape = sign), data=filtered_joined_clip_merged, color='red') +
  scale_shape_manual('FDR:', values = c(8, 16, 17, 1)) +
  xlab('log2FC (24°C/6°C) transcriptome') +
  ylab('log2FC (24°C/6°C) proteome') +
  theme_light() +
  guides(shape = guide_legend(override.aes = list(size=4))) +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15))

ggsave(file.path(dir,
                 'Eve_wo_EveB24_4',
                 paste0('transcr_proteome_logfc_', species,
                             '_FDR_24Cvs6C.png')), 
       scale = 1.1, width = 12, height = 8)

#### Boxplots of sign changed transcripts and corresponded proteins ####

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
                    species, '/count') # 3h

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
                    '6C_rep1', '6C_rep2', '6C_rep3') #,'6C_rep4') # Eve wo EveB24_4
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
combined_data$geneSymbol <- sub(', partial|LOW QUALITY PROTEIN:', '', combined_data$geneSymbol)
combined_data$geneSymbol <- sub('isoform X[0-9]', '', combined_data$geneSymbol)
var_width <- 23
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
         '> 0.05',
  ifelse(combined_data$sign == "< 0.05 (transcriptome)" & combined_data$method == 'MS/MS',
         '> 0.05', '< 0.05')))

combined_data$method <- factor(combined_data$method, 
                               levels = c('RNAseq', 'MS/MS'),
                               labels = c('Transcript', 'Protein'))
combined_data$condition <- factor(combined_data$condition, 
                                  levels = c('6C', '24C'),
                                  labels = c('6 °C (control)', '24.6 °C (heat stress)'))
combined_data_up <- subset(combined_data, log2FoldChange > 0) # up for transcriptomics
# or #
combined_data_up <- subset(combined_data, log2FoldChange > 0)
combined_data_up <- subset(combined_data_up, sign == '< 0.05 (transcriptome)' | (sign == '< 0.05 (both)' & logFC < 0))
write.table(combined_data_up, 
            file = file.path(dir,
                             'Eve_wo_EveB24_4',
                             paste(species, '_AllupDEtranscripts_joinedWithProteins.csv')))
combined_data_down <- subset(combined_data, log2FoldChange < 0) # down for transcriptomics
# or #
combined_data_down <- subset(combined_data, log2FoldChange < 0)
combined_data_down <- subset(combined_data_down, sign == '< 0.05 (transcriptome)' | (sign == '< 0.05 (both)' & logFC > 0))
write.table(combined_data_down, 
            file = file.path(dir, 
                             'Eve_wo_EveB24_4',
            paste(species, '_AlldownDEtranscripts_joinedWithProteins.csv')))

f <- function(x) {
  sapply(strsplit(x, '|', fixed=T), `[`, 2)
}

fdr_threshold <- 0.01
lfc_threshold <- 3.5
toplot <- subset(combined_data_up, pvalue_tran < fdr_threshold)
toplot <- subset(combined_data_down, pvalue_tran < fdr_threshold)
toplot <- subset(toplot, abs(log2FoldChange) > lfc_threshold)
toplot <- subset(toplot, geneSymbol != '*')
#toplot$sign2plot <- sub('0.05', fdr_threshold, toplot$sign2plot)

toplot <- subset(toplot, !grepl('uncharacterized protein|hypothetical', geneSymbol))

ggplot(toplot, aes(x = method, y = values, color = condition)) +
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
        strip.text = element_text(size = 8, margin = margin(2,2,2,2)),
        legend.position="bottom", 
        legend.margin = margin(0,0,1,0))

dir_to_save <- '~/labeglo2/proteome_transcr_comparision/slice_min_pvalue/'
dir_to_save <- '~/labeglo2/proteome_transcr_comparision/Eve_wo_EveB24_4/'
ggsave(file.path(dir_to_save, paste0(tolower(species), 
                                     '_DEtranscriptsWITHnotchangedProteins_down_fdr',
                                     fdr_threshold, '_lfc', lfc_threshold,
                                     '.png')),
       scale = 0.8,
       width = 9, height = 11.5)

#ggsave(file.path(dir_to_save, paste0(tolower(species), 
#                                     '_DEproteinsWITHtranscripts_up.pdf')),
#       scale = 0.8,
#       width = 12, height = 7)
# Ecy: width = 12.5, height = 7
# Gla: width = 9.5, height = 5
# Eve: width = 12.5, height = 7

### GO term analysis

library(topGO)
library(readr)
# upload go annotation file
gene_to_go <- read_delim('~/labeglo2/MS_results/annotation/go_gaf/combined.gaf',
                         '\t', comment = '!', col_names = F, 
                         col_types = cols(.default = "c"))
gene_to_go_bonus <- read_delim('~/labeglo2/MS_results/annotation/go_gaf/combined_bonus.csv',
                               '\t') # output from goterms_enrichment.py
# create named list with genes and corresponded go terms:
gene_go <- aggregate(X5 ~ X3, gene_to_go, function(x) list(unique(x)))
gene_go$X4 <- toupper(gene_go$X3)
gene_go_ <- gene_go$X5
names(gene_go_) <- gene_go$X4

# create named the same list but for the proteins with alternative names (bonus from goterms_enrichment.py) 
gene_go_bonus <- aggregate(go_term ~ gene_name, gene_to_go_bonus, 
                           function(x) list(unique(x)))
gene_go_bonus$up_gene_name <- toupper(gene_go_bonus$gene_name)
gene_go_bonus_ <- gene_go_bonus$go_term
names(gene_go_bonus_) <- gene_go_bonus$up_gene_name

# create the dataframe with proteins and information to them:
lfc_threshold <- 0
pvalue_threshold <- 0.01
unique_proteins <- joined_clip_merged$geneSymbol
joined_clip_merged$to_describe <- 
  ifelse(joined_clip_merged$log2FoldChange > lfc_threshold & joined_clip_merged$sign == '< 0.05 (transcriptome)' & joined_clip_merged$pvalue_tran <= pvalue_threshold,
         1, 0) # basically, make a boolean list (transcipts needed to describe vs all others)
joined_clip_merged$to_describe <- 
  ifelse(joined_clip_merged$log2FoldChange < lfc_threshold & joined_clip_merged$sign == '< 0.05 (transcriptome)' & joined_clip_merged$pvalue_tran <= pvalue_threshold,
         1, 0)
### take wanted statistics 
proteins_list <- joined_clip_merged$to_describe

### name vector with statistics by the corresponded protein:
names(proteins_list) <- toupper(unique_proteins)
proteins_list2 <- proteins_list[names(proteins_list) != '*'] # avoid not annotated proteins
proteins_list3 <- proteins_list2[names(proteins_list2) %in% gene_go$X4] # take only proteins with known go terms

# enrich the analysis by the observations with alternative names:
proteins_not_in_main_list <- proteins_list2[!names(proteins_list2) %in% gene_go$X4]
gene_go_from_bonus <- gene_go_bonus_[names(gene_go_bonus_) %in%
                                       names(proteins_not_in_main_list)]

# update go term list
gene_go_all_ <- c(gene_go_, gene_go_from_bonus)

# update the list with proteins in analysis
proteins_list4 <- proteins_list2[names(proteins_list2) %in% names(gene_go_all_)]
proteins_list4 <- factor(proteins_list4)
### run topgo
#thresh <- 0.01
ontology <- 'BP'
GOdata <- new("topGOdata", description = "Simple session", 
              ontology = ontology,
              allGenes = proteins_list4, 
              #geneSel = function(x) x < 0.01, # pMM, pGS
              #geneSel = function(x) abs(x) > threshold, # MM, GS
              #geneSel = function(x) x < thresh, # MM
              nodeSize = 8,
              annot = annFUN.gene2GO, 
              gene2GO = gene_go_all_)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", 
                   topNodes = 70,
                   numChar = 100)
allGO <- genesInTerm(GOdata)

go_table_w_genes <- NULL

for (row in 1:nrow(allRes)) {
  current_row <- allRes[row,]
  current_go <- as.character(current_row[1])
  wanted_ <- proteins_list4[names(proteins_list4) %in% allGO[[current_go]]] 
  genes_names <- paste0(names(wanted_[wanted_ == 1]), collapse = '', sep = ';')
  names(genes_names) <- 'significant_genes'
  go_table_w_genes <- rbind(go_table_w_genes, 
                            as.data.frame(c(current_row, genes_names)))
}

write.csv(go_table_w_genes,
          file = file.path(dir_to_save, 
                           paste0(species, 
                                  '_TOPGO_GOterms_TRANSCRIPTS_withNOchangeInProteins_',
                                  'DOWN', 
                                  '_lfc', lfc_threshold,
                                  '_FDR', pvalue_threshold,
                                  '.csv')))


head(annot[match(rownames(counts), annot$contig),]$annot)
counts_df <- data.frame(counts)
counts_df$annotation <- annot[match(rownames(counts), annot$contig),]$annot
transcr$annot <- annot[match(transcr$contig, annot$contig),]$annot
