#!/bin/bash

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load Trimmomatic/0.39-Java-1.8

##########################################################################################

# Don't forget to edit the files directory

Prefix=$(sed -n ${SLURM_ARRAY_TASK_ID}p Data/Dataset.txt | cut -f1) # To run as an array

# FileName R1
RNA_f1=$(grep "$Prefix.*_R1_" ${path_data}List_RNAseq_files_fullpath.txt)
# Number of FileName R1 (== number of sequencing run)
wc_l_f1=$(grep "$Prefix.*_R1_" ${path_data}List_RNAseq_files_fullpath.txt | wc -l)

# FileName R2
RNA_f2=$(grep "$Prefix.*_R2_" ${path_data}List_RNAseq_files_fullpath.txt)
# Number of FileName R2 (== number of sequencing run)
wc_l_f2=$(grep "$Prefix.*_R2_" ${path_data}List_RNAseq_files_fullpath.txt | wc -l)



if [ $wc_l_f1 -eq 1 ]; then

        echo 'One sequencing run'
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 4 \
        /${RNA_f1} /${RNA_f2} \
        ${Prefix}_trimmed_paired_1.fastq.gz ${Prefix}_trimmed_unpaired_1.fastq.gz ${Prefix}_trimmed_paired_2.fastq.gz ${Prefix}_trimmed_unpaired_2.fastq.gz \
        ILLUMINACLIP:${adaptors}:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:80

else
        echo 'Two sequencing runs'

        RNA_f1_tmp1=$(echo $RNA_f1 | cut -d' ' -f1)
        RNA_f1_tmp2=$(echo $RNA_f1 | cut -d' ' -f2)

        RNA_f2_tmp1=$(echo $RNA_f2 | cut -d' ' -f1)
        RNA_f2_tmp2=$(echo $RNA_f2 | cut -d' ' -f2)


        RNA_f1_comb=${Prefix}_R1_combined.fastq
        RNA_f2_comb=${Prefix}_R2_combined.fastq

        if [ ! -f $RNA_f1_comb -a ! -f $RNA_f2_comb ]; then

                less ${RNA_f1_tmp1} > $RNA_f1_comb
                less ${RNA_f1_tmp2} >> $RNA_f1_comb

                less ${RNA_f2_tmp1} > $RNA_f2_comb
                less ${RNA_f2_tmp2} >> $RNA_f2_comb

        fi

        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 4 \
        $RNA_f1_comb $RNA_f2_comb \
        ${Prefix}_trimmed_paired_1.fastq.gz ${Prefix}_trimmed_unpaired_1.fastq.gz ${Prefix}_trimmed_paired_2.fastq.gz ${Prefix}_trimmed_unpaired_2.fastq.gz \
        ILLUMINACLIP:${adaptors}:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:80

fi

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
