#### To look at proteins with all zeros in one condition ######

library(tidyverse) 
library(limma) 
library(edgeR)
library(ggfortify) 

species <- 'Gla'

### 1. Upload metafile
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity corrected', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  meta <- meta[grepl(species_name, meta$condition),]
  return(meta)
}
path2meta <-
  'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/Metadata_Proteus.tsv'

meta <- meta_upload(path2meta, species)

# if you want to take all 6C controls as one control group:
meta$condition <- ifelse(grepl('CK1|CK2|VrK1|VrK2|LK1|LK2', meta$sample), 
                         '6C', '24C')

### 2. Data uploading (proteinGroups file from MaxQuant)
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_eve/' # 24h
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_ecy/' # 24h
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_gla/' # 24h
proteinGroups_file <- 'proteinGroups_wo_cont_more2pept.txt' # take file with proteinGroups with 2 or more peptides quantified
dat_init <- read.csv(file.path(dir, proteinGroups_file), sep = '\t', header = T, 
                     check.names = F) 
pep_annot <- read.csv(file.path(dir, 
                                paste0('annot_protein_groups_', tolower(species),'.csv')), 
                      sep = '\t',
                      header = T)

select_data <- function(meta_data, proteinGroups_data){
  dat_ <- proteinGroups_data[,c('Protein IDs', meta_data$measure)]
  rownames(dat_) <- dat_$`Protein IDs`
  dat_ <- dat_[-1]
  colnames(dat_) <- meta_data$sample
  dat_[dat_ == 0] <- NA
  return(dat_)
}
dat <- select_data(meta, dat_init)

### 3. Get rid of the outliers and samples, which you do not want to analyze
dat <- dat[, !grepl('VrK1_390_3|VrK1_390_4|pool', colnames(dat))]
dat <- dat[, !grepl('pool', colnames(dat))]
meta <- subset(meta, sample != 'VrK1_390_3' & sample != 'VrK1_390_4' &
                 !grepl('pool', sample))
meta <- subset(meta, !grepl('pool', sample))

condition <- as.factor(meta$condition)
ixs <- order(condition)
dat_ <- dat[, ixs]
condition <- condition[ixs]

for (level in levels(condition)) {
  col_ixs <- condition == level
  dat_ <- dat_[rowSums(is.na(dat_[,col_ixs])) == sum(col_ixs),]
}

View(dat_)
# For E. verrucosus: 15 proteins with NAs in all samples -> exclude from analysis
# For E. cyaneus: 11 proteins with NAs in all samples -> exclude from analysis
# For G. lacustris: 9 proteins with NAs in all samples -> exclude from analysis

# It turned out that we do not have proteins with all zeros in only one condition



