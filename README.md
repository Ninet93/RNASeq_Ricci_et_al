# Original article: Visual opsin gene expression dynamics in the adaptive radiation of cichlid fishes of Lake Tanganyika


## Authors Information

* Virginie Ricci (virginie.ricci@unibas.ch)
    * Department of Environmental Sciences, Zoological Institute, University of Basel, Basel, Switzerland

* Fabrizia Ronco
    * Natural History Museum, University of Oslo, Oslo, Norway

* Nicolas Boileau
    * Department of Environmental Sciences, Zoological Institute, University of Basel, Basel, Switzerland

* Walter Salzburger (walter.salzburger@unibas.ch)
    * Department of Environmental Sciences, Zoological Institute, University of Basel, Basel, Switzerland


## How to cite this article

Ricci, V., Ronco, F., Boileau, N., & Salzburger, W. (in revision). Visual opsin gene expression dynamics in the adaptive radiation of cichlid fishes of Lake Tanganyika



## Overview

In this study, we sequenced 753 new retinal transcriptomes of 112 cichlid species from African Lake Tanganyika to (i) identify adaptive changes in the expression of rod and cone visual opsin genes, (ii) define and reconstruct the evolution of visual palettes on the basis of cone opsin expression levels, and (iii) examine rod and cone opsin expression levels in relation to macro-habitat, diet, and relative eye size.
The raw reads are available from NCBI under the BioProject accession no. PRJNA913112 (https://www. ncbi.nlm.nih.gov/bioproject/)


## Methodological information

We generated a retinal read count dataset of 753 cichlid individuals using _Nile tilapia_ as reference sequence.
We examined gene expression patterns in the retina using principal component analysis (PCA).
We identified opsin expression profiles and characterised visual palettes.
We investigated ecological and morphological correlates of visual opsin gene expression using linear regression (lm), phylogenetic generalized least squares (pGLS) and phylogenetic partial least square regression (pPLS)
Please, visit also the GitHub repository for the scripts used in this study (https://github.com/Ninet93/RNASeq_Ricci_et_al)



## Scripts (GitHub repository)

* 00_Parse_GTF_biotypes.sh: Genes information of _Nile tilapia_ RefSeq (RefSeq accession GCF_001858045.2)

* 01_Trimmomatic.sh: Trimmomatic

* 02_FastQC_MultiQC.sh: FastQC and MultiQC

* 03_STAR_genome_index.sh: STAR genome indexing

* 04_STAR_mapping.sh: STAR mapping for RNA-seq data

* 04_STAR_mapping_RH2As.sh: STAR mapping of unmapped reads to retrieve read pairs that map to both RH2Aa and RH2Ab

* 05_MappedReadPairsToRH2As.sh: Get only read pairs that map to both RH2Aa and RH2Ab

* 06_1_HTSeqCount.sh: HTSeq-count (count reads in genomic features)

* 06_2_HTSeqCount_RH2As.sh: HTSeq-count of read pairs that map to both RH2Aa and RH2Ab
 
* 07_HTSeqCount_results.R: Remove lowly expressed genes, select protein-coding RNAs and lncRNAs, convert read count to TPM, create DESEq2 objects

* 08_HTSeqCount_results_WeightedSpeciesMean_opsins.R: Calculate weighted species mean TPM values

* 09_DESeq2_results.R: Principal Component Analysis (PCA) and phylospace on PC

* 010_VisualPalettes.R: Visual palettes and ancestral states reconstructions

* 011_Correlations.R: Linear regression (lm), phylogenetic generalized least squares (pGLS) and two-block partial least squares regression (PLS) analyses
  
* 012_Cocorrelations.R: Co-correlation of visual opsin gene expression levels




## Data (GitHub repository)

### Inputs

* Orenil_opsins.txt: _Nile tilapia_ opsins information
	* 3 columns: Visual opsin gene, GeneName, Color
	* 8 rows
	
* RNAseq_SpeciesTree.tre: Species tree (see Ronco et al. 2021) pruned to the taxa used in this study

### Outputs

* GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt: _Nile tilapia_ genome annotation summary tables generated using FeatureCount (FeatureCount is not used in this study)
	* 6 columns: GeneID, Chromosome, Start, End, Strand, Length
	* 41947 rows

* GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt: _Nile tilapia_ genes information tables
	* 3 columns: Scaffold, GeneID, BioType
	* 42623 rows

* HTSeqCount_ALL_exons_Individuals_withRH2As.txt: Raw HTSeq-count table with the correction of RH2Aa/RH2Ab proportion (row=individuals, column=genes)
	* 41946 columns: GeneID
	* 754 rows

### Extra Figures

* HTSeqCount_exons_PCA_Individual_proteins_lncRNAs.pdf: Individual PCA of protein-coding RNAs and lncRNAs

* HTSeqCount_exons_PCA_Individual_opsins.pdf: Individual PCA of opsins

* HTSeqCount_RodCones_PE_Count_exons_PerTribe.pdf: Individual barplots of opsins (HTSeq-count)

* HTSeqCount_RodCones_PE_TPM_exons_PerTribe.pdf: Individual barplots of opsins (HTSeq-count converted in TPM values)


#### This README was created by V. Ricci (2023)