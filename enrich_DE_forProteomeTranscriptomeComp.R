### To add Hemocyanin subunit 1 contig from transcriptome of Ecy: 
# At first, run signDEprot_with_RNAseqData.R until inner_join with proteins and transcripts data
# The similarity between two contigs (Gla and Ecy) ~95%
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_4259_length_3812_cov_87.630880_g2442_i0'] <- 
  'NODE_13910_length_2189_cov_1646.189252_g8989_i0'
# Then proceed with signDEprot_with_RNAseqData.R

### 14-3-3zeta protein for Gla
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_31737_length_1118_cov_2319.587465_g6963_i3'] <- 
  'NODE_46667_length_850_cov_2110.021223_g6133_i3'
# The similarity between two contigs (Eve and Gla) ~90%, actually they form the same protein group

### Superoxide dismutase copper/zinc binding domain for Eve
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_15269_length_1950_cov_481.873751_g8862_i1'] <- 
  'NODE_36202_length_988_cov_849.140575_g19385_i2'

### MAN2C1 for Eve
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_30038_length_1235_cov_5.642496_g20401_i0'] <- 
  'NODE_7788_length_2991_cov_9.625765_g4914_i0'

### hemocyanin subunit 1 for Eve:
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_4259_length_3812_cov_87.630880_g2442_i0'] <- 
  'NODE_12577_length_2250_cov_4404.478873_g8154_i0'

### arylsulfatase I-like for Eve:
# only 72% of similarity with Eve sequence
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_36737_length_1037_cov_13.795547_g25241_i0'] <- 
  'NODE_4552_length_3904_cov_219.468223_g2822_i0'

### eukaryotic translation initiation factor 2 subunit 1-like for Eve:
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_9042_length_2655_cov_101.478511_g5359_i0'] <- 
  'NODE_16839_length_1860_cov_282.252899_g11015_i0'

### Pgm for Eve:
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_12575_length_2326_cov_157.593764_g8101_i0'] <- 
  'NODE_6500_length_3303_cov_60.894284_g4084_i0'

### Peptidase C1A propeptide for Eve:
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_13197_length_2139_cov_171.834928_g8058_i0'] <- 
  'NODE_13617_length_2140_cov_215.235294_g8884_i0'


###### Hsps
### HSPA1A-inducible for Ecy
# the similarity is 98.27% with the original Eve variant;
# change straight after target_info table creation
target_info$protein[which(target_info$protein ==
        'NODE_3822_length_4220_cov_2696.563654_g2376_i0.p1_Eve')[2]] <- 
  'NODE_5378_length_3701_cov_2676.348850_g3383_i0.p1_Ecy'
# change only the second one (for Ecy)
