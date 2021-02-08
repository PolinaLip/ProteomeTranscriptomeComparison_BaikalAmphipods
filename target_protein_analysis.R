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

to_plot <- subset(combined_data, orthogroup == 'OG0000002')
 
ggplot(to_plot, aes(x = condition, y = intensities)) +
  geom_boxplot(aes(color = species), outlier.alpha = 0) +
  geom_jitter(aes(color = species)) +
  scale_color_manual('Species', values = c('#7570b3', '#1b9e77', '#d95f02')) +
  facet_wrap(~species + protein_group + annotation + hsp70_type, scales = 'free') +
  theme_bw()



