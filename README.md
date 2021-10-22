## ProteomeTranscriptomeComparison_BaikalAmphipods
Here are scripts for the comparison of the proteome and transcriptome of the Baikal amphipods exposed to 24C (for 24 hours, heat shock).  

### prot_transcr_lfc.R: 
1. relates each protein group to their transcripts; 
2. performs correlation test;
3. plots correlation plots between LFC from transcriptome DE analysis and LFC from proteome DE analysis.

### zerosCond_heatshock.R: 
I checked proteins for cases when there are all NAs in one condition. We do not have such proteins but we have about ten proteins with all NAs in all condition. I removed theprot_transcr_lfc.R: m from the analysis.
