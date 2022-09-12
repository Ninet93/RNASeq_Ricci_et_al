#!/bin/bash

#SBATCH --job-name=Parse_GTF
#SBATCH --mail-user=virginie.ricci@unibas.ch
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --output=/scicore/home/salzburg/fityxu76/VisualSystem/Pipeline/Logs/RNA_seq/Parse_GTF_%A.out
#SBATCH --error=/scicore/home/salzburg/fityxu76/VisualSystem/Pipeline/Logs/RNA_seq/Parse_GTF_%A.err
#SBATCH --workdir=/scicore/home/salzburg/fityxu76/VisualSystem/Pipeline/PE_RNAseq/Athimed/

DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of job' $DATE

module purge
module load Python/3.5.2-goolf-1.7.20

#################################################################################################################################
# combSJs.py coded by Geff - SciCORE
#################################################################################################################################

path_Onil='/scicore/home/salzburg/GROUP/NCBI_Orenil_GCF_001858045_2/ncbi-genomes-2020-01-13/'
gtf_f='GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gtf'

path_VR='/scicore/home/salzburg/GROUP/NCBI_Orenil_GCF_001858045_2/ncbi-genomes-2020-01-13/VRicci'
#gtf_f=test.gtf

python Parse_GTF_proteins_lncRNAs_test.py $gtf_f $path_Onil $path_VR

#python Parse_GTF.py $gtf_f $path_VR $path_VR


#################################################################################################################################
#################################################################################################################################

DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of job' $DATE
