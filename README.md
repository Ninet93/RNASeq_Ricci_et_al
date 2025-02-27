# Visual opsin gene expression dynamics in the adaptive radiation of cichlid fishes of Lake Tanganyika
Virginie Ricci, Fabrizia Ronco, Nicolas Boileau & Walter Salzburger (2023)

## Scripts

* [`00_Parse_GTF_biotypes.sh`](Scripts/00_Parse_GTF_biotypes.sh): Genes information of _Nile tilapia_ RefSeq (RefSeq accession GCF_001858045.2)
  * Python and Parse_GTF_biotypes.py

* [`01_Trimmomatic.sh`](Scripts/01_Trimmomatic.sh): Trimmomatic
  * Trimmomatic

* [`02_FastQC_MultiQC.sh`](Scripts/02_FastQC_MultiQC.sh): FastQC and MultiQC
  * FastQC and MultiQC

* [`03_STAR_genome_index.sh`](Scripts/03_STAR_genome_index.sh): STAR genome indexing
  * STAR

* [`04_STAR_mapping.sh`](Scripts/04_STAR_mapping.sh): STAR mapping for RNA-seq data
  * STAR and SAMtools

* [`04_STAR_mapping_RH2As.sh`](Scripts/04_STAR_mapping_RH2As.sh): STAR mapping of unmapped reads to retrieve read pairs that map to both RH2Aa and RH2Ab
  * STAR and SAMtools

* [`05_MappedReadPairsToRH2As.sh`](Scripts/05_MappedReadPairsToRH2As.sh): Get only read pairs that map to both RH2Aa and RH2Ab
  * SAMtools

* [`06_1_HTSeqCount.sh`](Scripts/06_HTSeqCount.sh): HTSeq-count (count reads in genomic features)
  * HTSeq and SAMtools

* [`06_2_HTSeqCount_RH2As.sh`](Scripts/06_HTSeqCount_RH2As.sh): HTSeq-count of read pairs that map to both RH2Aa and RH2Ab
  * HTSeq and SAMtools
 
* [`07_HTSeqCount_results.R`](Scripts/HTSeqCount_results.R): Remove lowly expressed genes, select protein-coding RNAs and lncRNAs, convert read count to TPM, create DESEq2 objects
  * R, DESeq2

* [`08_HTSeqCount_results_WeightedSpeciesMean_opsins.R`](Scripts/08_HTSeqCount_results_WeightedSpeciesMean_opsins.R): Calculate weighted species mean TPM values
  * R

* [`09_DESeq2_results.R`](Scripts/09_DESeq2_results.R): Principal Component Analysis (PCA) and phylospace on PCA
  * R, DESeq2

* [`010_VisualPalettes.R`](Scripts/010_VisualPalettes.R): Visual palettes and ancestral states reconstructions
  * R

* [`011_Correlations.R`](Scripts/011_Correlations.R): Linear regression (lm), phylogenetic generalized least squares (pGLS) and two-block partial least squares regression (PLS) analyses
  * R
  
* [`012_Cocorrelations.R`](Scripts/012_Cocorrelations.R): Expression correlation of visual opsin gene expression levels
  * R





## Data

### Inputs

* [`Orenil_opsins.txt`](Data/Orenil_opsins.txt): _Nile tilapia_ opsins information

* [`RNAseq_SpeciesTree.tre`](Data/RNAseq_SpeciesTree.tre): Species tree (see Ronco et al. 2021) pruned to the taxa used in this study

### Outputs

* [`GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt`](Data/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt): _Nile tilapia_ genome annotation summary tables generated using FeatureCount (FeatureCount is not used in this study)

* [`GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt`](Data/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt): _Nile tilapia_ genes information tables

* [`HTSeqCount_ALL_exons_Individuals_withRH2As.txt`](Data/HTSeqCount_ALL_exons_Individuals_withRH2As.txt): Raw HTSeq-count table with the correction of RH2Aa/RH2Ab proportion (row=individuals, column=genes)

### Extra Figures

* [`HTSeqCount_exons_PCA_Individual_proteins_lncRNAs.pdf`](Data/Extra_figures/HTSeqCount_exons_PCA_Individual_proteins_lncRNAs.pdf): Individual PCA of protein-coding RNAs and lncRNAs

* [`HTSeqCount_exons_PCA_Individual_opsins.pdf`](Data/Extra_figures/HTSeqCount_exons_PCA_Individual_opsins.pdf): Individual PCA of opsins

* [`HTSeqCount_RodCones_PE_Count_exons_PerTribe.pdf`](Data/Extra_figures/HTSeqCount_RodCones_PE_Count_exons_PerTribe.pdf): Individual barplots of opsins (HTSeq-count)

* [`HTSeqCount_RodCones_PE_TPM_exons_PerTribe.pdf`](Data/Extra_figures/HTSeqCount_RodCones_PE_TPM_exons_PerTribe.pdf): Individual barplots of opsins (HTSeq-count converted in TPM values)

* [`Opsins_median_readcoverage_CDS_GW_PerTribe.pdf`](Data/Extra_figures/Opsins_median_readcoverage_CDS_GW_PerTribe.pdf): Coverage of mapped read estimates per cone opsin (see Ricci et al. 2022 for RH1)
