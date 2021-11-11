## ProteomeTranscriptomeComparison_BaikalAmphipods
Here are scripts for the comparison between the proteome and transcriptome of two endemic Baikal amphipods (E. verrucosus and E. cyaneus) and a Holarctic species (G. lacustris) exposed to 24C (for 24 hours, heat shock experiment). The here studied species are non-model organisms with no full-assembled genome available.

### prot_transcr_lfc.R: 
1. relates each protein group to their transcripts; 
2. performs correlation test;
3. plots correlation plots between LFC from transcriptome DE analysis and LFC from proteome DE analysis.

### zerosCond_heatshock.R: 
I checked proteins for cases when there are all NAs in one condition. We do not have such proteins but we have about ten proteins with all NAs in all condition. I removed them from the analysis.

### signDEprot_with_RNAseqData.R: 
To plot boxplots with MS/MS and RNAseq data for all differentialy abundant proteins.

### fig_cor.R:
To plot fig. 1 - Correlation between log2 of fold changes of transcriptome and proteome at 24.6 °C compared to 6 °C (24h proteome with 24h and 3h transcriptome).
