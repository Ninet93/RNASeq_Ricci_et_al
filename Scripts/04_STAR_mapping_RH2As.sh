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

# Don't forget to edit the files directory

Prefix=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f1) # To run as an array

# RefSeq fasta and genome annotation files (GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna from NCBI)
FASTA=GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna
GTF=GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gtf

RNA_f1=${Prefix}_Unmapped.out.mate1 # unmapped reads from first round of STAR
RNA_f2=${Prefix}_Unmapped.out.mate2 # unmapped reads from first round of STAR

# Retrieve reads that map exactly to RH2Aa and RH2Ab (percent identity of exons between 94.81% and 99.40%)

STAR --runMode alignReads --genomeDir $PATH --outFileNamePrefix ${Prefix}_Unmapped_ --readFilesIn $RNA_f1 $RNA_f2 --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
 --outFilterMultimapNmax 2 --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 \
 --outWigType bedGraph --outWigStrand Stranded --outWigNorm None --outSAMattrRGline ID:${Prefix} PL:ILLUMINA LB:NovaSeq_6000


##########################################################################################

module purge
module load SAMtools/1.10-foss-2018b

##########################################################################################

samtools index -b ${Prefix}_Unmapped_Aligned.sortedByCoord.out.bam

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE