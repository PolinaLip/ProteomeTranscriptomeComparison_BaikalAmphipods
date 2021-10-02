##### to add ORF type to the data and select only complete ORFs and then proceed with 
##### de_separately.R with data_irs 
pr <- data.frame(pr_names = rownames(data_irs),
                 pr_gr = 1:length(rownames(data_irs)))

pr_sep <-  separate(pr, pr_names, ';', into=paste0('key', 1:30), 
                    fill='right')
pr_long <- pivot_longer(pr_sep, starts_with('key'),
                        names_to = NULL, values_to = 'protein') %>%
           filter(!is.na(protein))

dir <- '~/labeglo2/proteome_transcr_comparision/'
orf_type <- read.table(file = file.path(dir, 'all_species_orf_type.tsv'))
colnames(orf_type) <- c('protein_name', 'orf', 'length')

pr_long$orftype <- 
  orf_type[match(pr_long$protein, orf_type$protein_name),]$orf

pr_long_filt <- pr_long %>%
  group_by(pr_gr) %>%
  filter(orftype == 'complete') 

inx <- unique(pr_long_filt$pr_gr)

#data_irs_ <- data_irs
data_irs <- data_irs[inx,]
