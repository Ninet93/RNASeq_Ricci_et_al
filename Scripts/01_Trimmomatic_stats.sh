#!/bin/bash

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE
##########################################################################################

STAT=Trimmomatic_stats.txt
touch $STAT
rm $STAT
echo $'ID\tNbSeqRuns\tReadPairs\tSurvivingBoth\tSurvivingBoth_Percent\tSurvivingF\tSurvivingF_Percent\tSurvivingR\tSurvivingR_Percent\tDropped\tDropped_Percent' > $STAT

# 756 individuals
for nb in {1..756}; do
        OUT=Trimmomatic_PE_${ID_job}_${nb}.out
        ERR=Trimmomatic_PE_${ID_job}_${nb}.err

        ID=$(sed -n 2p $OUT) # ID sample
        nb_seq=$(sed -n 4p $OUT) # number of sequencing run

        readpairs=$(sed -n 6p $ERR | cut  -d' ' -f4) # number of read pairs

        survivingreadpairs=$(sed -n 6p $ERR | cut  -d' ' -f7) # number of surviving pairs after trimming
        percent_survivingreadpairs=$(sed -n 6p $ERR | cut  -d' ' -f8 | sed 's/(//' | sed 's/)//')

        survivingF=$(sed -n 6p $ERR | cut  -d' ' -f12) # number of surviving forward reads after trimming
        percent_survivingF=$(sed -n 6p $ERR | cut  -d' ' -f13 | sed 's/(//' | sed 's/)//')

        survivingR=$(sed -n 6p $ERR | cut  -d' ' -f17) # number of surviving reverse reads after trimming
        percent_survivingR=$(sed -n 6p $ERR | cut  -d' ' -f18 | sed 's/(//' | sed 's/)//')

        dropped=$(sed -n 6p $ERR | cut  -d' ' -f20) # number of dropped reads after trimming
        percent_dropped=$(sed -n 6p $ERR | cut  -d' ' -f21 | sed 's/(//' | sed 's/)//')


        echo -e $ID'\t'$nb_seq'\t'$readpairs'\t'$survivingreadpairs'\t'$percent_survivingreadpairs'\t'$survivingF'\t'$percent_survivingF'\t'$survivingR'\t'$percent_survivingR'\t'$dropped'\t'$percent_dropped >> $STAT

done

##########################################################################################
DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE