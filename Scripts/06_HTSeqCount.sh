#!/bin/bash

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load HTSeq/0.11.2-foss-2018b-Python-3.6.6

##########################################################################################

# Don't forget to edit the files directory

Prefix=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f1) # To run as an array

GTF=GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf # manually corrected (GeneID for tRNAs, etc...)

htseq-count -f bam -r pos -m union -s reverse -t exon ${Prefix}_Aligned_PE.sortedByCoord.out.bam $GTF > ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt_tmp

tail -5 ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt_tmp > ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt_summary

head -n -5 ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt_tmp > ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt_0

sed 1d ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt_0 > ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt

rm ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt_tmp ${Prefix}_Aligned_PE.sortedByCoord.out.count.txt_0


##########################################################################################

module purge
module load SAMtools/1.10-foss-2018b

##########################################################################################

samtools flagstat ${Prefix}_Aligned_PE.sortedByCoord.out.bam > ${Prefix}_Aligned_PE.sortedByCoord.out.bam.stats
 
##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE