#!/bin/bash

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load FastQC/0.11.8-Java-1.8

##########################################################################################

# Don't forget to edit the files directory

Prefix=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f1) # To run as an array


fastqc -o $PATH ${Prefix}_trimmed_paired_1.fastq.gz ${Prefix}_trimmed_paired_2.fastq.gz
fastqc -o $PATH ${Prefix}_trimmed_unpaired_1.fastq.gz ${Prefix}_trimmed_unpaired_2.fastq.gz

##########################################################################################

module purge
module load MultiQC/1.11-foss-2018b-Python-3.6.6

##########################################################################################

multiqc -d *_paired_*_fastqc.zip -o $PATH # Run once when fastqc jobs are done

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE