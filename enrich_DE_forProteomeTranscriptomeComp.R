### To add Hemocyanin protein 1 contig from transcriptome of Ecy: 
# At first, run signDEprot_with_RNAseqData.R until inner_join with proteins and transcripts data
# The similarity between two contigs (Gla and Ecy) ~95%
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_4259_length_3812_cov_87.630880_g2442_i0'] <- 
  'NODE_13910_length_2189_cov_1646.189252_g8989_i0'
# Then proceed in signDEprot_with_RNAseqData.R

### 14-3-3zeta protein for Gla
proteins_long$protein_clip[proteins_long$protein_clip == 
                             'NODE_31737_length_1118_cov_2319.587465_g6963_i3'] <- 
  'NODE_46667_length_850_cov_2110.021223_g6133_i3'
# The similarity between two contigs (Eve and Gla) ~90%, actually they form the same protein group

