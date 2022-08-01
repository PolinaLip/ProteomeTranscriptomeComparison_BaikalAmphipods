### to check the correlation between proteome and transcriptome (HS) using 
### TPM for transcripts and "intensities" values for proteins (rather than "intensity count")
library(tidyverse)
library(ggplot2)

trans_dir <- '~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/Eve/counts/'
transcripts <- data.frame()
for (filename in Sys.glob(file.path(trans_dir, '*/quant.sf'))) {
  sample <- dirname(strsplit(filename, '//')[[1]][2])
  transcripts <- rbind(
      transcripts,
                         
      read.csv(filename, sep = '\t') %>%
        select(c('Name', 'TPM')) %>%
        add_column(Sample = sample)
    )
}

#transcripts <- pivot_wider(transcripts, names_from = 'Sample', values_from = 'TPM') %>%
#  as.data.frame() # if you want scale the data by the row
#rownames(transcripts) <- transcripts$Name # if you want scale the data by the row
#transcripts <- transcripts[-1] # if you want scale the data by the row
#transcripts_colnames <- colnames(transcripts) # if you want scale the data by the row
#transcripts_scaled <- apply(transcripts, 1, scale) # if you want scale the data by the row
#transcripts_scaled <- t(transcripts_scaled) # if you want scale the data by the row
#colnames(transcripts_scaled) <- transcripts_colnames # if you want scale the data by the row
#transcripts_scaled <- data.frame(transcripts_scaled) # if you want scale the data by the row
#transcripts <- transcripts_scaled %>%
#  rownames_to_column(var = 'Name') %>%
#  pivot_longer(!Name,
#               names_to = 'Sample', values_to = 'TPM') # if you want scale the data by the row

transcripts <- mutate(transcripts,
                      Condition = ifelse(grepl('_24C_', Sample), '24C', '6C'))

prot_dir <- '~/labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_eve'
intensities <- read.csv(file.path(prot_dir, 'intensities_after_slNorm_eve.csv'),
                        sep = '\t')
#rownames(intensities) <- intensities$protein_group # if you want scale the data by the row
#intensities <- intensities[-1] # if you want scale the data by the row
#names_of_columns <- colnames(intensities) # if you want scale the data by the row
#intensities <- apply(intensities, 1, scale) # if you want scale the data by the row
#intensity_scaled <- t(intensities) # if you want scale the data by the row
#intensity_scaled <- data.frame(intensity_scaled) # if you want scale the data by the row
#colnames(intensity_scaled) <- names_of_columns # if you want scale the data by the row

intensities <- pivot_longer(!protein_group,
               names_to = 'Sample', values_to = 'Intensity')

#intensities <- intensity_scaled %>% # if you want scale the data by the row
#  rownames_to_column(var = 'protein_group') %>% # if you want scale the data by the row
#  pivot_longer(!protein_group, # if you want scale the data by the row
#               names_to = 'Sample', values_to = 'Intensity') # if you want scale the data by the row

intensities <- mutate(intensities,
      Condition = ifelse(grepl('K', Sample), '6C', '24C'),
      protein_group = gsub('\\.p[0-9]+_(Eve|Gla|Ecy)', '', protein_group))

trans_aggr <- aggregate(TPM ~ Name + Condition, transcripts, mean)
intens_aggr <- aggregate(Intensity ~ protein_group + Condition, intensities, mean) %>%
  mutate(Name = protein_group)

# shared_names <- sort(intersect(intens_aggr$protein_group, trans_aggr$Name))
# trans_aggr2 <- trans_aggr[match(shared_names, trans_aggr$Name),]
# intens_aggr2 <- intens_aggr[match(shared_names, intens_aggr$protein_group),]
joined <- inner_join(trans_aggr, intens_aggr, by = c('Name', 'Condition'))
joined_control <- filter(joined, Condition == '6C' & TPM != 0)
joined_treatment <- filter(joined, Condition == '24C' & TPM != 0)

with(joined_control, cor.test(log(TPM), log(Intensity)))
with(joined_treatment, cor.test(log(TPM), log(Intensity)))
with(joined_control, cor.test(TPM, Intensity))
with(joined_treatment, cor.test(TPM, Intensity))

ggplot(joined_control) +
  geom_point(aes(TPM, Intensity)) +
  theme_bw()
