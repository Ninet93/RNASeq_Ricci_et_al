# Visual opsin expression profile in the adaptive radiation of cichlid fishes of Lake Tanganyika
Virginie Ricci, Fabrizia Ronco, Nicolas Boileau & Walter Salzburger (2022)

## Scripts

* [`00_Parse_GTF_biotypes.sh`](Scripts/00_Parse_GTF_biotypes.sh): Genes information of Nile tilapia RefSeq (RefSeq accession GCF_001858045.2)
  * Python and Parse_GTF_biotypes.py

* [`01_Trimmomatic.sh`](Scripts/01_Trimmomatic.sh): Trimmomatic
  * Trimmomatic

* [`01stats_Trimmomatic.sh`](Scripts/01stats_Trimmomatic.sh): Trimmomatic summary statistics

* [`02_FastQC_MultiQC.sh`](Scripts/02_FastQC_MultiQC.sh): FastQC and MultiQC

* [`03_STAR_genome_index.sh`](Scripts/03_STAR_genome_index.sh): STAR genome indexing

* [`04_STAR_mapping.sh`](Scripts/04_STAR_mapping.sh): STAR mapping for RNA-seq data

* [`04_STAR_mapping_RH2As.sh`](Scripts/04_STAR_mapping_RH2As.sh): STAR mapping of unmapped reads to retrieve read pairs that map to both RH2Aa and RH2Ab

* [`05_MappedReadPairsToRH2As.sh`](Scripts/05_MappedReadPairsToRH2As.sh): Get only read pairs that map to both RH2Aa and RH2Ab

* [`06_HTSeqCount.sh`](Scripts/06_HTSeqCount.sh): HTSeq-count (count reads in genomic features)

* [`06_HTSeqCount_RH2As.sh`](Scripts/06_HTSeqCount_RH2As.sh): HTSeq-count of read pairs that map to both RH2Aa and RH2Ab
 





## Data

### Inputs



### Outputs

* [`GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt`](Data/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt): Nile tilapia genome annotation summary tables generated using FeatureCount (FeatureCount is not used in this study)

* [`GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt`](Data/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt): Nile tilapia genes information tables
