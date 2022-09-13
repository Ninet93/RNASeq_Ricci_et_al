#!/bin/bash

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

module purge
module load SAMtools/1.7-goolf-1.7.20

##########################################################################################

# Don't forget to edit the files directory

# Get reads mapped to RH2Aa and RH2Ab
# RH2Aa
samtools view -Sb ${Prefix}_Unmapped_Aligned.sortedByCoord.out.bam "NC_031970.2:16274618-16285731" > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.sorted.bam_tmp
samtools view -h $${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.sorted.bam_tmp | awk 'substr($0,1,1)=="@" || ($4>= 16274618 && $4<= 16285731)' | samtools view -Sb > $${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.sorted.bam
rm $${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.sorted.bam_tmp

# RH2Ab
samtools view -Sb ${Prefix}_Unmapped_Aligned.sortedByCoord.out.bam "NC_031970.2:16286976-16289027" > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.sorted.bam_tmp
samtools view -h ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.sorted.bam_tmp | awk 'substr($0,1,1)=="@" || ($4>= 16286976 && $4<= 16289027)' | samtools view -Sb > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.sorted.bam
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.sorted.bam_tmp

for GENE in RH2Aa RH2Ab; do
        samtools index ${Prefix}_Unmapped_Aligned.sortedByCoord.out.${GENE}.sorted.bam
done

##########################################################################################

# Get reads ID for RH2Aa and RH2Ab
samtools view ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.sorted.bam | cut -f1 | sort | uniq > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID.txt

samtools view ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.sorted.bam | cut -f1 | sort | uniq > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID.txt

# Get common reads ID
comm -12 ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID.txt ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID.txt > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.common_ID.txt

# Get number of common reads ID
nb_common=$(wc -l ${Prefix}_Unmapped_Aligned.sortedByCoord.out.common_ID.txt | cut -d' ' -f1)

# Get BAM lines of common reads ID
touch ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID.bam
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID.bam
touch ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID.bam
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID.bam
for x in $(seq 1 $nb_common); do
        ID=$(sed -n ${x}p ${Prefix}_Unmapped_Aligned.sortedByCoord.out.common_ID.txt)
        samtools view ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.sorted.bam | grep $ID >> ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID.bam
        samtools view ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.sorted.bam | grep $ID >> ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID.bam
done


# Get alignment of mulimapped reads - potentially in RH2Aa and RH2Ab only
grep 'NH:i:2' ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID.bam > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.bam
grep 'NH:i:2' ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID.bam > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.bam

# get read pairs only
cut -d$'\t' -f1 ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.bam | sort | uniq -c | grep '2 A*' > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.count
cut -d$'\t' -f1 ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.bam | sort | uniq -c | grep '2 A*' > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.count

sed -i 's/      2 A/A/g' ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.count
sed -i 's/      2 A/A/g' ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.count

sort ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.count > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.count.sorted
sort ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.count > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.count.sorted

# Get common read pairs only
comm -12 ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.count.sorted ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.count.sorted > RH2Aa_RH2Ab_common.txt

# Get number of common read pairs only
nb_l=$(wc -l RH2Aa_RH2Ab_common.txt | cut -d' ' -f1)

# Get BAM lines of common read pairs only
touch ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.bam_tmp
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.bam_tmp
touch ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.bam_tmp
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.bam_tmp
for x in $(seq 1 $nb_l); do
        id=$(sed -n ${x}p RH2Aa_RH2Ab_common.txt)

        grep $id ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.bam >> ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.bam_tmp
        grep $id ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.bam >> ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.bam_tmp

done

rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID.txt
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID.txt
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.common_ID.txt
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID.bam
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID.bam
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.bam
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.bam
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.count
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.count
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.count.sorted
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.count.sorted
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa_ID_firstalign.bam_tmp
rm ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab_ID_firstalign.bam_tmp


##########################################################################################

module purge
module load SAMtools/1.15-GCC-10.3.0

##########################################################################################

# Get BAM lines of common read pairs only in BAM format

samtools view -hN RH2Aa_RH2Ab_common.txt ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.sorted.bam | samtools view -Sb > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Aa.common.bam
samtools view -hN RH2Aa_RH2Ab_common.txt ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.sorted.bam | samtools view -Sb > ${Prefix}_Unmapped_Aligned.sortedByCoord.out.RH2Ab.common.bam

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE