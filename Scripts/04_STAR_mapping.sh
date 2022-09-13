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

STAR --runMode alignReads --genomeDir $PATH --outFileNamePrefix ${Prefix}_ --readFilesIn $RNA_f1 $RNA_f2 --readFilesCommand zcat --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
 --outFilterMultimapNmax 1 --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 \
 --outWigType bedGraph --outWigStrand Stranded --outWigNorm None --outSAMattrRGline ID:${Prefix} PL:ILLUMINA LB:NovaSeq_6000


##########################################################################################

module purge
module load SAMtools/1.10-foss-2018b

##########################################################################################

samtools index -b ${Prefix}_Aligned.sortedByCoord.out.bam

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE