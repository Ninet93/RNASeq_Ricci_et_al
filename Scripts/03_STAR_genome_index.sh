#!/bin/bash

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load STAR/2.7.3a-foss-2018b

##########################################################################################

# RefSeq fasta and genome annotation files (GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna from NCBI)
FASTA=GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna
GTF=GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gtf

STAR --runMode genomeGenerate --genomeDir $PATH --sjdbGTFfile $GTF --sjdbOverhang 99  --genomeFastaFiles $FASTA --runThreadN 8

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE