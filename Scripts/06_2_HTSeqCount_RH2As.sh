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

samtools merge ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.common.bam ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.common.bam ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.common.bam

htseq-count -f bam -r pos -m union --nonunique all -a 0 -o ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.common.count.sam -s reverse -t exon ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.common.bam $GTF > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt_tmp

tail -5 ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt_tmp > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt_summary

head -n -5 ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt_tmp > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt_0

sed 1d ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt_0 > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt

rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt_tmp ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.count.txt_0


##########################################################################################

module purge
module load SAMtools/1.10-foss-2018b

##########################################################################################

samtools flagstat ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.common.bam > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_RH2Ab.common.bam.stats
 
##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE