# counts data preparation
meta_upload <- function(path_to_file, species_name) {
  meta <- read.csv(file = path_to_file, sep = '\t')
  meta$measure <- sub('intensity', 'intensity count', meta$measure)
  meta$measure <- paste(meta$measure, meta$experiment)
  meta <- meta[grepl(species_name, meta$condition),]
  return(meta)
}
path2meta <-
  'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/Metadata_Proteus.tsv'
meta <- meta_upload(path2meta, 'Eve')

dir <- 'labeglo2/MS_results/390/withDBfromRNAspades/wIMBR2/protein_groups_eve/'
proteinGroups_file <- 'proteinGroups_wo_cont_more2pept.txt' # take file with proteinGroups with 2 or more peptides quantified
dat_init <- read.csv(file.path(dir, proteinGroups_file), sep = '\t', header = T, 
                     check.names = F) 

select_data <- function(meta_data, proteinGroups_data){
  dat <- proteinGroups_data[,c('Protein IDs', meta_data$measure)]
  rownames(dat) <- dat$`Protein IDs`
  dat <- dat[-1]
  colnames(dat) <- meta$sample
  dat[dat == 0] <- NA
  return(dat)
}
dat_count <- select_data(meta, dat_init)
dat_count <- dat_count[,!grepl('VrK1_390_3|VrK1_390_4', colnames(dat_count))]
dat_count <- dat_count[,!grepl('Eve_pool', colnames(dat_count))]
meta <- subset(meta, sample != 'VrK1_390_3' & sample != 'VrK1_390_4' &
                 !grepl('Eve_pool', sample)
)

proteins_intensities_no_na <- na.omit(dat_count)
proteins_intensities_no_na <- as.data.frame(proteins_intensities_no_na)
data_raw_count <- proteins_intensities_no_na
boxplot(log2(data_raw_count), col = c(rep(c('red','blue'), each = 9), rep('green', 4)), 
        notch = TRUE, main = 'RAW data: Exp1 (red), Exp2 (blue), Exp3 (green)',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_raw_count), col = c(rep(c('red','blue'), each = 9), rep('green', 4)), 
              main = 'Raw data', legend=F) # Eve, resolution -> 800/553

exp1_raw_count <- data_raw_count[c(1:8)]
exp2_raw_count <- data_raw_count[c(9:16)]
exp3_raw_count <- data_raw_count[c(17:19)]

target_count <- mean(c(colSums(exp1_raw_count, na.rm = T), 
                 colSums(exp2_raw_count, na.rm = T),
                 colSums(exp3_raw_count, na.rm = T)), na.rm = T) # Eve 
norm_facs_count <- target_count / colSums(exp1_raw_count, na.rm = T)
exp1_sl_count <- sweep(exp1_raw_count, 2, norm_facs_count, FUN = "*")
norm_facs_count <- target_count / colSums(exp2_raw_count, na.rm = T)
exp2_sl_count <- sweep(exp2_raw_count, 2, norm_facs_count, FUN = "*")
norm_facs_count <- target_count / colSums(exp3_raw_count, na.rm = T)
exp3_sl_count <- sweep(exp3_raw_count, 2, norm_facs_count, FUN = "*")

data_sl_count <- cbind(exp1_sl_count, exp2_sl_count, exp3_sl_count)

boxplot(log2(data_sl_count), col = c(rep(c('red', 'blue'), each = 9), rep('green', 4)), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nExp1 (red), Exp2 (blue), Exp3 (green)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity') # Eve, resolution -> 1400/690
plotDensities(log2(data_sl_count), col = c(rep(c('red', 'blue'), each = 9), rep('green', 4)),
              main = "SL normalization", legend = F) # Eve, resolution -> 800/553

plotMDS(log2(data_sl_count), col = c(rep(c('red', 'blue'), each = 9), rep('green', 4)), 
        main = "SL clusters group by TMT experiment") # Eve
data_sl_t_count <- na.omit(data_sl_count) %>% t %>% as.data.frame
pca_res <- prcomp(data_sl_t_count)
meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
autoplot(pca_res, data=meta2, colour='condition')

irs <- tibble(rowSums(exp1_sl_count), 
              rowSums(exp2_sl_count), rowSums(exp3_sl_count)) 
colnames(irs) <- c("sum1", "sum2", "sum3") 
irs$average <- apply(irs, 1, function(x) exp(mean(log(x), na.rm = T)))

irs$fac1 <- irs$average / irs$sum1 * ncol(exp1_sl_count) # Ecy, Eve, Gla
irs$fac2 <- irs$average / irs$sum2 * ncol(exp2_sl_count) # Ecy, Eve, Gla
irs$fac3 <- irs$average / irs$sum3 * ncol(exp3_sl_count) # Eve, Gla

data_irs_count <- exp1_sl_count * irs$fac1
data_irs_count <- cbind(data_irs_count, exp2_sl_count * irs$fac2)
data_irs_count <- cbind(data_irs_count, exp3_sl_count * irs$fac3)

boxplot(log2(data_irs_count), col = c(rep(c('red', 'blue'), each = 9), rep('green', 4)), 
        main = "Internal Reference Scaling (IRS) normalized data: \nExp1 (red), Exp2 (blue), Exp3 (green)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity', notch = TRUE) # Eve, resolution 1400/690

col_vec <- ifelse(grepl('pool', colnames(data_irs_count)), 'green',
                  ifelse(grepl('VrK1', colnames(data_irs_count)), 'orange',
                         ifelse(grepl('Vr1', colnames(data_irs_count)), 'blue', 'red')))

plotMDS(log2(data_irs_count), col = col_vec, main = "IRS clusters do not group by TMTsets anymore")

data_irs_t_count <- na.omit(data_irs_count) %>% t %>% as.data.frame
pca_res <- prcomp(data_irs_t_count)
meta2 <- meta
meta2$experiment <- as.factor(meta2$experiment)
autoplot(pca_res, data=meta2, colour='condition')

write.csv(data_irs_count, 
          file = file.path('labeglo2/proteome_transcr_comparision/norm_counts_Eve_HS_protGroups.csv'))

