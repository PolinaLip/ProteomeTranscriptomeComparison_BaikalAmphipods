## ProteomeTranscriptomeComparison_BaikalAmphipods
Here are scripts for the comparison between the proteome and transcriptome of two endemic Baikal amphipods (E. verrucosus and E. cyaneus) and a Holarctic species (G. lacustris) exposed to 24C (for 24 hours, heat shock experiment). The here studied species are non-model organisms with no full-assembled genome available.

### proteinGroup_annotation.py:
To combine output files from DIAMOND and EggNOG annotations and create an annotation file with all necessary information for protein groups. This annotation file is used in the data analysis scripts in R (DE_analysis_[species_name].R, signDEprot_with_RNAseqData.R, etc.).

### combine_annotations.py:
To combine output files from DIAMOND and EggNOG annotations and create an annotation file with all necessary information for transcripts. This annotation file is used in the data analysis scripts in R (prot_transcr_conc.R, prot_transcr_lfc.R, signDEprot_with_RNAseqData.R, SIGNtranscr_nonSIGNprot.R)

### DE_analysis_[species_name].R:
For differential expression analysis of the transcriptomes of different species using "counts" tables and DESeq2. 

### de_separately.R:
For differential abundance (DA) analysis of proteins performed by EdgeR. We performed DA analysis separately for each number of missing observations in rows because the edgeR package forbids data with missing values. The missing values were most probably arisen due to selectivity of data-dependent acquisition and batch-effect caused by batches of TMT-labeling.

### prot_transcr_lfc.R: 
1. relates each protein group to their transcripts; 
2. performs correlation test;
3. plots correlation plots between LFC from transcriptome DE analysis and LFC from proteome DE analysis.

### zerosCond_heatshock.R: 
I checked proteins for cases when there are all NAs in one condition. We do not have such proteins but we have about ten proteins with all NAs in all condition. I removed them from the analysis.

### signDEprot_with_RNAseqData.R: 
To plot boxplots with MS/MS and RNAseq data for all differentialy abundant proteins.
## enrich_DE_forProteomeTranscriptomeComp.R:
For several DA proteins, the protein group names contained only contig names with assignment to other species (e.g. in the case of E. verrucosus proteome, there were protein groups contains contigs with _Gla (contig came from transcriptome assembly of G. lacustris)). I found the closest sequences for these protein groups in the transcriptome assembly of the targeted species; and I used these transcipts for the analysis and the plotting. 

### fig_cor.R:
To plot fig. 1 - Correlation between log2 of fold changes of transcriptome and proteome at 24.6 °C compared to 6 °C (24h proteome with 24h and 3h transcriptome).
