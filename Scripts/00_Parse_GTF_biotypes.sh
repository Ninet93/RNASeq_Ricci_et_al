#!/bin/bash

##########################################################################################
### Virginie Ricci, 2021
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load Python/3.5.2-goolf-1.7.20

##########################################################################################

# Don't forget to edit the files directory

path_Onil='/NCBI_Orenil_GCF_001858045_2/ncbi-genomes-2020-01-13/'
gtf_f='GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gtf'

python Parse_GTF_biotypes.py $gtf_f $path_Onil

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE
