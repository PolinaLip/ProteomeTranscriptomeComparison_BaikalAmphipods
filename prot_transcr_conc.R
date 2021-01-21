library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(rstatix)

### Intensity values preparation:
intensities <- read.csv('labeglo2/proteome_transcr_comparision/norm_counts_Ecy_HS_protGroups_withNA.csv')
rownames(intensities) <- intensities$X
intensities <- intensities[-1]
intensities <- log2(intensities)
meta <- 
  read.csv(file = 
             'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/Metadata_Proteus.tsv',
           sep = '\t')
meta <- meta[grepl('Eve', meta$condition),]
meta <- subset(meta, sample != 'VrK1_390_3' & sample != 'VrK1_390_4')

intensities <- intensities[,grepl('VrK2|Vr1|CK2|C1|LK2|L1', 
                                  colnames(intensities))]

pep_annot <- 
  read.csv('labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_ecy/annot_protein_groups_ecy.csv', 
           sep = '\t',
           header = T)

prot2annot <- match(rownames(intensities), pep_annot$protein_group)
sum(is.na(prot2annot))
intensities$annotation <- pep_annot$annotation[prot2annot]

intens_long <- pivot_longer(intensities %>% rownames_to_column('protein_group'),
                            starts_with(c('C', 'V', 'L')), names_to = 'sample', 
                            values_to = 'intensity')
intens_long$condition <- meta$condition[match(intens_long$sample, meta$sample)]
intens_long$condition <- ifelse(grepl('K', intens_long$sample), '6C', '24C')

intens_long <- as.data.frame(intens_long)
intens_long_mean <- aggregate(intensity ~ protein_group + condition + annotation, 
                              data = intens_long, median, na.rm = T)

### Counts preparation:
annot <- read.csv('labeglo2/proteome_transcr_comparision/contigs_whole_annot_Ecy.csv', 
                  sep = '\t')
dir_transcr <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/Ecy/counts' # specify path to your samples

files_names <- c('Eve_6C_rep1', 'Eve_6C_rep2', 'Eve_6C_rep3', 'Eve_6C_rep4', 
                 'Eve_24C_rep1', 'Eve_24C_rep2', 'Eve_24C_rep3', 'Eve_24C_rep4') # specify the names of folders with quant.sf data, Eve

files_names <- c('Ecy_6C_rep1', 'Ecy_6C_rep2', 'Ecy_6C_rep3', 'Ecy_6C_rep4', 
                 'Ecy_24C_rep1', 'Ecy_24C_rep2', 'Ecy_24C_rep3', 'Ecy_24C_rep4') # specify the names of folders with quant.sf data, Ecy

files_names <- c('Gla_6C_rep1', 'Gla_6C_rep2', 'Gla_6C_rep3', 'Gla_6C_rep4', 
                 'Gla_24C_rep1', 'Gla_24C_rep2', 'Gla_24C_rep3') # specify the names of folders with quant.sf data, Gla

countFiles <- list.files(dir_transcr, full.names = T)
countFiles
counts <- lapply(countFiles, function(countsFile) {
  f <- file.path(countsFile, 'quant.sf')
  read.table(f, sep="\t", header=1, row.names = 1, 
             stringsAsFactors = F, comment.char = "#")
})
counts <- lapply(counts, function(countsTable) countsTable[, 3, drop=F])

counts <- do.call(cbind, counts)
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3', '24C_rep4',
                  '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4')
colnames(counts) <- sample_names

contig2annot <- match(rownames(counts), annot$contig)
counts$transcr_annot <- annot$annot[contig2annot]

counts_long <- pivot_longer(counts %>% rownames_to_column('contig'), 
                            starts_with(c('2','6')), names_to = 'transc_sample',
                            values_to = 'counts')

counts_long$condition <- ifelse(grepl('24C', counts_long$transc_sample), 
                                '24C', '6C')
counts_long_mean <- aggregate(counts ~ contig + condition + transcr_annot, 
                              data = counts_long, median, na.rm = T)

### 24C and 6C correlation

intens_long_mean_24C <- subset(intens_long_mean, condition == 'Gla_24C')
intens_long_mean_24C <- intens_long_mean_24C[-2]
intens_long_mean_6C <- subset(intens_long_mean, condition == 'Gla_6C_after')
intens_long_mean_6C <- subset(intens_long_mean, condition == '6C')
intens_long_mean_6C <- intens_long_mean_6C[-2]

# choose condition:
intens_to_analise <- intens_long_mean_24C
intens_to_analise <- intens_long_mean_6C

intens_to_analise$prot_group_id <- 1:nrow(intens_to_analise)
proteins_sep <- separate(intens_to_analise, protein_group, ';', 
                         into=paste0('key', 1:30), 
                         fill='right')
proteins_long <- pivot_longer(proteins_sep, starts_with('key'),
                              names_to = NULL, values_to = 'protein') %>%
  filter(!is.na(protein))

# make a column with protein names without everything after .p (including .p) 
proteins_long$protein_clip <- sub('\\.p[0-9]+_...$', '', proteins_long$protein)

# if there are two proteins with the same clipped name take the protein with 
# maximum intensity value:
proteins_long <- group_by(proteins_long, protein_clip) %>%
  slice_max(intensity, n = 1, with_ties = F)

counts_long_mean_24C <- subset(counts_long_mean, condition == '24C')
counts_long_mean_24C <- counts_long_mean_24C[-2]
counts_long_mean_6C <- subset(counts_long_mean, condition == '6C')
counts_long_mean_6C <- counts_long_mean_6C[-2]

# choose condition:
counts_to_analyse <- counts_long_mean_24C
counts_to_analyse <- counts_long_mean_6C

joined <- inner_join(counts_to_analyse, proteins_long, 
                     by = c('contig' = 'protein_clip'))

joined_merged <- joined %>% 
  group_by(prot_group_id) %>%
  mutate(max_count = max(counts)) %>% 
  mutate(transcr_groups = paste0(contig, collapse = ';')) %>%
  mutate(proteome_groups = paste0(protein, collapse = ';')) %>%
  ungroup() %>%
  dplyr::select(!c(contig, counts, protein, transcr_annot)) %>%
  unique()

ggplot(joined, aes(log2(counts), intensity)) +
  geom_point(alpha=.5, color='gray70') +
  geom_smooth() +
#  geom_point(data=filtered_joined_clip, color='red') +
#  geom_text_repel(aes(label= geneSymbol),
#                  data = filtered_joined_clip) +
#  xlab('log2FC (24°C/6°C) transcriptome') +
#  ylab('log2FC (24°C/6°C) proteome') +
  theme_light() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))
# this plot will give duplicated values, below is without duplicates
cor.test(joined_merged$max_count, exp(joined_merged$intensity))

tmp <- filter(proteins_long, grepl('Serpin', annotation))
ggplot(joined_merged, aes(log2(max_count), intensity)) +
  geom_point(aes(color=grepl(paste0(tmp$protein_clip, collapse='|'), transcr_groups)), alpha=.5)
tmp2 <- filter(joined_merged,
               grepl(paste0(tmp$protein_clip, collapse='|'), transcr_groups))

cor_test_res <- cor.test(joined_merged$max_count, exp(joined_merged$intensity))

joined_merged <- subset(joined_merged, max_count != 0)

ggplot(joined_merged, aes(log2(max_count), intensity)) +
  geom_point(alpha=.5, color='gray70') +
  geom_smooth(method='lm') +
  #  geom_point(data=filtered_joined_clip, color='red') +
  #  geom_text_repel(aes(label= geneSymbol),
  #                  data = filtered_joined_clip) +
  annotate(geom='text', x = -4.5, y = 12, hjust = 0,
           label = paste0('r2 = ', round(cor_test_res$estimate, 4), '\n',
                          'p-value ', p_format(cor_test_res$p.value,
                                                 accuracy = 0.000001))) +
  xlab('transcipt abundance at 6°C (log2 count)') +
  ylab('protein abundance at 6°C (log2 count)') +
  theme_light() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave('labeglo2/proteome_transcr_comparision/prot_trans_counts_Ecy_6Callcontrols_withNAs_lm.png',
       scale = 1.5)

filter(joined, log2(counts) > 10 & intensity < 0.1)
tmp <- filter(intens_long, grepl('NODE_1149_length_6503_cov_6439', protein_group))

# use Reporter intensity count from proteinGroup.txt
# if use intensities with all proteins (even with zero inside) -> we can see some positive correlation between proteome and transcriptome
#ggsave('labeglo2/proteome_transcr_comparision/prot_trans_conc_ProteinsWithZerosAlso.png')

# ! use only 6C from parallel control ! -> and the same applies to lfc comparision of proteome and transcriptome

