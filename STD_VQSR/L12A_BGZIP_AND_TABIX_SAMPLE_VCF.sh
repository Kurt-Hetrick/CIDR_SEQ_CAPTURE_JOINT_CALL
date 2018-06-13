#########################################
##### ---qsub parameter settings--- #####
#########################################

###########################################################
### --THESE ARE OVERWRITTEN IN THE PIPELINE DURING QSUB ###
###########################################################

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge to submit any of these queue when available
#$ -q prod.q,rnd.q

# tell sge that you are in the users current working directory
#$ -cwd

# tell sge to export the users environment variables
#$ -V

# tell sge to submit at this priority setting
#$ -p -1000

# tell sge to output both stderr and stdout to the same file
#$ -j y

#######################################
##### END QSUB PARAMETER SETTINGS #####
#######################################

# export all variables, useful to find out what compute node the program was executed on
set

# INPUT VARIABLES

TABIX_DIR=$1

CORE_PATH=$2
PROJECT_SAMPLE=$3
SM_TAG=$4

START_BGZIP_AND_TABIX=`date '+%s'`

CMD1=$TABIX_DIR'/bgzip -c '$CORE_PATH'/'$PROJECT_SAMPLE'/VCF/RELEASE/FILTERED_ON_BAIT/'$SM_TAG'_MS_OnBait.vcf'
CMD1=$CMD1' >| '$CORE_PATH'/'$PROJECT_SAMPLE'/VCF/RELEASE/FILTERED_ON_BAIT/TABIX/'$SM_TAG'_MS_OnBait.vcf.gz'

CMD2=$TABIX_DIR'/tabix -f -p vcf '$CORE_PATH'/'$PROJECT_SAMPLE'/VCF/RELEASE/FILTERED_ON_BAIT/TABIX/'$SM_TAG'_MS_OnBait.vcf.gz'

END_BGZIP_AND_TABIX=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_SAMPLE",L01,BGZIP_AND_TABIX,"$HOSTNAME","$START_BGZIP_AND_TABIX","$END_BGZIP_AND_TABIX \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".WALL.CLOCK.TIMES.csv"

echo $CMD1 >> $CORE_PATH/$PROJECT_SAMPLE/command_lines.txt
echo >> $CORE_PATH/$PROJECT_SAMPLE/command_lines.txt
echo $CMD1 | bash

echo $CMD2 >> $CORE_PATH/$PROJECT_SAMPLE/command_lines.txt
echo >> $CORE_PATH/$PROJECT_SAMPLE/command_lines.txt
echo $CMD2 | bash

