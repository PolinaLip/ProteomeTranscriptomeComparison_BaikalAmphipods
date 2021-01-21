###########################################
#### Proteome counts data preparation #####
###########################################

library(limma)
library(dplyr)
library(tidyverse)
library(ggfortify)

### 0. Choose species
species <- 'Eve'

### 1. Upload metafile
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity count', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  meta <- meta[grepl(species_name, meta$condition),]
  return(meta)
}
path2meta <-
  'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/Metadata_Proteus.tsv'

meta <- meta_upload(path2meta, species)

### 2. Upload proteinGroups file 
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_eve/'
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_ecy/'
dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_gla/'

proteinGroups_file <- 'proteinGroups_wo_cont_more2pept.txt' # take file with proteinGroups with 2 or more peptides quantified
dat_init <- read.csv(file.path(dir, proteinGroups_file), sep = '\t', header = T, 
                     check.names = F) 

### 3. Choose columns with counts data
select_data <- function(meta_data, proteinGroups_data){
  dat <- proteinGroups_data[,c('Protein IDs', meta_data$measure)]
  rownames(dat) <- dat$`Protein IDs`
  dat <- dat[-1]
  colnames(dat) <- meta$sample
  dat[dat == 0] <- NA
  return(dat)
}

dat_count <- select_data(meta, dat_init)

### 4. Remove samples which are outliers (should be determined in advance) 
###    or just unwanted one
dat_count <- dat_count[,!grepl('VrK1_390_3|VrK1_390_4', colnames(dat_count))]
dat_count <- dat_count[,!grepl('pool', colnames(dat_count))]
meta <- subset(meta, sample != 'VrK1_390_3' & sample != 'VrK1_390_4' &
                 !grepl('pool', sample))

#proteins_intensities_no_na <- na.omit(dat_count)
#proteins_intensities_no_na <- as.data.frame(proteins_intensities_no_na)
#data_raw_count <- proteins_intensities_no_na

### 5. To get rid of observations with all NAs in one of the conditions
cond <- as.factor(meta$condition)
ixs <- order(cond)
dat_count <- dat_count[, ixs]
cond <- cond[ixs]

for (level in levels(cond)) {
  col_ixs <- cond == level
  dat_count <- dat_count[rowSums(is.na(dat_count[,col_ixs])) < sum(col_ixs),]
}

data_raw_count <- dat_count

### 6. Sample loading and IRS normalization

sl_irs_normalization <- function(counts, metafile, species) {
  if (species == 'Eve'){
    exp1_raw <- counts[,subset(metafile, experiment == 6)$sample]
    exp2_raw <- counts[,subset(metafile, experiment == 7)$sample]
    exp3_raw <- counts[,subset(metafile, experiment == 8)$sample]
  }
  else if (species == 'Ecy') {
    exp1_raw <- counts[,subset(metafile, experiment == 4)$sample]
    exp2_raw <- counts[,subset(metafile, experiment == 5)$sample]
  }
  else if (species == 'Gla') {
    exp1_raw <- counts[,subset(metafile, experiment == 1)$sample]
    exp2_raw <- counts[,subset(metafile, experiment == 2)$sample]
    exp3_raw <- counts[,subset(metafile, experiment == 3)$sample]
  }
  if (species == 'Eve' | species == 'Gla') {
    target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T), 
                     colSums(exp3_raw, na.rm = T)),
                   na.rm = T)
    norm_facs <- target / colSums(exp1_raw, na.rm = T)
    exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
    norm_facs <- target / colSums(exp2_raw, na.rm = T)
    exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
    norm_facs <- target / colSums(exp3_raw, na.rm = T)
    exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")
    
    data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl)
    
    irs <- tibble(rowSums(exp1_sl), 
                  rowSums(exp2_sl), rowSums(exp3_sl)) 
    colnames(irs) <- c("sum1", "sum2", "sum3") 
    irs$average <- apply(irs, 1, function(x) exp(mean(log(x), na.rm = T)))
    
    irs$fac1 <- irs$average / irs$sum1 * ncol(exp1_sl) 
    irs$fac2 <- irs$average / irs$sum2 * ncol(exp2_sl) 
    irs$fac3 <- irs$average / irs$sum3 * ncol(exp3_sl) 
    
    data_irs_count <- exp1_sl * irs$fac1
    data_irs_count <- cbind(data_irs_count, exp2_sl * irs$fac2)
    data_irs_count <- cbind(data_irs_count, exp3_sl * irs$fac3)
  }
  
  else if (species == 'Ecy') {
    target <- mean(c(colSums(exp1_raw, na.rm = T), colSums(exp2_raw, na.rm = T)),
                   na.rm = T)
    norm_facs <- target / colSums(exp1_raw, na.rm = T)
    exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
    norm_facs <- target / colSums(exp2_raw, na.rm = T)
    exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
    
    data_sl <- cbind(exp1_sl, exp2_sl)
    
    irs <- tibble(rowSums(exp1_sl), rowSums(exp2_sl)) 
    colnames(irs) <- c("sum1", "sum2") 
    irs$average <- apply(irs, 1, function(x) exp(mean(log(x), na.rm = T)))
    
    irs$fac1 <- irs$average / irs$sum1 * ncol(exp1_sl) 
    irs$fac2 <- irs$average / irs$sum2 * ncol(exp2_sl) 
    
    data_irs_count <- exp1_sl * irs$fac1
    data_irs_count <- cbind(data_irs_count, exp2_sl * irs$fac2)
  }
  return(data_irs_count)
}

data_sl_irs_count <- sl_irs_normalization(data_raw_count, 
                                          meta, 
                                          species = species)

boxplot(log2(data_sl_irs_count), notch = TRUE,
        xlab = 'TMT Sample', ylab = 'log2 of counts') # Eve, resolution -> 1400/690
plotDensities(log2(data_sl_irs_count), 
              main = "SL normalization", legend = F) # Eve, resolution -> 800/553

plotMDS(log2(data_sl_irs_count), col = c(rep(c('red', 'blue'), each = 9), 
                                     rep('green', 4)), 
        main = "SL clusters group by TMT experiment") # Eve

data_sl_irs_count_t <- na.omit(data_sl_irs_count) %>% t %>% as.data.frame
pca_res <- prcomp(data_sl_irs_count_t)
meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
autoplot(pca_res, data=meta2, colour='condition')

write.csv(data_sl_irs_count, 
          file = file.path('labeglo2/proteome_transcr_comparision/norm_counts_Eve_HS_protGroups_withNA.csv'))

# then -> prot_transcr_conc.R
