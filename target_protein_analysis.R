library(ggplot2)

target_info <- read.csv('~/labeglo2/MS_results/390/withDBfromRNAspades/hsps_orthologes/orthologues_hsps_sum.tsv',
                      sep = '\t') 
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

meta <- read.csv(file = 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/Metadata_Proteus.tsv',
                 sep = '\t')
meta$condition <- ifelse(grepl('CK1|CK2|VrK1|VrK2|LK1|LK2', meta$sample), 
                         '6C', '24C')

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

combined_data <- NULL
for (row in 1:nrow(target_info)) {
  if (target_info[row, 'species'] == 'Ecy'){
    combined_data <- 
      rbind(combined_data,
        connect_intens_and_info(target_info[row,], ecy_intensities))
  }  
  else if (target_info[row, 'species'] == 'Eve') {
    combined_data <- 
      rbind(combined_data,
            connect_intens_and_info(target_info[row,], eve_intensities))
  }
  else if (target_info[row, 'species'] == 'Gla') {
    combined_data <- 
      rbind(combined_data,
            connect_intens_and_info(target_info[row,], gla_intensities))
  }
}

colnames(combined_data)[7] <- 'intensities'
colnames(combined_data)[8] <- 'sample'
combined_data$condition <- ifelse(grepl('CK1|CK2|VrK1|VrK2|LK1|LK2', combined_data$sample), 
                    '6C', '24C')

to_plot <- subset(combined_data, orthogroup == 'OG0000022')
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
species <- 'Ecy'
species <- 'Eve'
species <- 'Gla'

dir <- paste0('~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/', 
              species, '/counts') # specify path to your samples
current_dir <- paste0('~/labeglo2/Transcriptomics/quantification/HS_exp_labeglo1/',
                      species)

countFiles <- list.files(dir, full.names = T)
countFiles
counts <- lapply(countFiles, function(countsFile) {
  f <- file.path(countsFile, 'quant.sf')
  read.table(f, sep="\t", header=1, row.names = 1, 
             stringsAsFactors = F, comment.char = "#")
})
counts <- lapply(counts, function(countsTable) countsTable[, 4, drop=F])

counts <- do.call(cbind, counts)
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3', '24C_rep4',
                  '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4') # Ecy, Eve
sample_names <- c('24C_rep1', '24C_rep2', '24C_rep3',
                  '6C_rep1', '6C_rep2', '6C_rep3', '6C_rep4') # Gla
colnames(counts) <- sample_names

counts_ecy <- counts
counts_eve <- counts
counts_gla <- counts

data_prep_counts <- function(data_counts){
  data_counts_scaled <- apply(data_counts, 1, scale)
  data_counts_scaled <- t(data_counts_scaled)
  colnames(data_counts_scaled) <- colnames(data_counts)
  return(data_counts_scaled)
}

counts_ecy <- data_prep_counts(counts_ecy)
counts_eve <- data_prep_counts(counts_eve)
counts_gla <- data_prep_counts(counts_gla)

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

colnames(combined_data_transc)[7] <- 'counts'
colnames(combined_data_transc)[8] <- 'sample'
combined_data_transc$condition <- ifelse(grepl('6C', combined_data_transc$sample), 
                                  '6C', '24C')

to_plot <- subset(combined_data_transc, orthogroup == 'OG0000000')
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
  facet_wrap(~all_labels, labeller = as_labeller(f), ncol = 4) +
  ylab('Scaled counts') +
  xlab('') +
  theme_bw() +
  theme(legend.text = element_text(face = 'italic', size = 11),
        axis.title.y = element_text(size = 12))
        #,
        #strip.text = element_text(size=7))

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



