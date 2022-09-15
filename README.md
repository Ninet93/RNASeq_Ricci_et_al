# Visual opsin expression profile in the adaptive radiation of cichlid fishes of Lake Tanganyika
Virginie Ricci, Fabrizia Ronco, Nicolas Boileau & Walter Salzburger (2022)

## Scripts

* [`00_Parse_GTF_biotypes.sh`](Scripts/00_Parse_GTF_biotypes.sh): Genes information of Nile tilapia RefSeq (RefSeq accession GCF_001858045.2)
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

* [`09_DESeq2_results.R`](Scripts/09_DESeq2_results.R): Principal Component Analysis (PCA)
  * R, DESeq2




## Data

### Inputs

* [`Orenil_opsins.txt`](Data/Orenil_opsins.txt): Nile tilapia opsin information

* [`RNAseq_SpeciesTree.tre`](Data/RNAseq_SpeciesTree.tre): Species tree (see Ronco et al. 2021)

### Outputs

* [`GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt`](Data/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt): Nile tilapia genome annotation summary tables generated using FeatureCount (FeatureCount is not used in this study)

* [`GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt`](Data/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt): Nile tilapia genes information tables
