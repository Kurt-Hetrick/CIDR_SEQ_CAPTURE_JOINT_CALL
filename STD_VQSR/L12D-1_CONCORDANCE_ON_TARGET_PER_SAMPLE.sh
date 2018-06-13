#########################################
##### ---qsub parameter settings--- #####
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

CIDR_SEQSUITE_JAVA_DIR=$1
CIDR_SEQSUITE_6_1_1_DIR=$2
VERACODE_CSV=$3

CORE_PATH=$4
PROJECT_SAMPLE=$5
SM_TAG=$6
TARGET_BED=$7

START_CONCORDANCE=`date '+%s'`

# decompress the filtered on target... 

zcat $CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_SNV.vcf.gz" \
$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG/$SM_TAG"_MS_OnTarget_SNV.vcf"

CMD=$CIDR_SEQSUITE_JAVA_DIR'/java -jar'
CMD=$CMD' '$CIDR_SEQSUITE_6_1_1_DIR'/CIDRSeqSuite.jar'
CMD=$CMD' -pipeline -concordance'
CMD=$CMD' '$CORE_PATH'/'$PROJECT_SAMPLE'/TEMP/'$SM_TAG
CMD=$CMD' '$CORE_PATH'/'$PROJECT_SAMPLE'/Pretesting/Final_Genotyping_Reports/'
CMD=$CMD' '$CORE_PATH'/'$PROJECT_SAMPLE'/TEMP/'$SM_TAG
CMD=$CMD' '$TARGET_BED
CMD=$CMD' '$VERACODE_CSV

echo $CMD >> $CORE_PATH/$PROJECT_SAMPLE/command_lines.txt
echo >> $CORE_PATH/$PROJECT_SAMPLE/command_lines.txt
echo $CMD | bash

END_CONCORDANCE=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_SAMPLE",M01,CONCORDANCE,"$HOSTNAME","$START_CONCORDANCE","$END_CONCORDANCE \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".WALL.CLOCK.TIMES.csv"

mv $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG/$SM_TAG"_concordance.csv" \
$CORE_PATH/$PROJECT_SAMPLE/REPORTS/CONCORDANCE_MS/$SM_TAG"_concordance.csv"

mv $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG/missing_data.csv \
$CORE_PATH/$PROJECT_SAMPLE/REPORTS/CONCORDANCE_MS/$SM_TAG"_missing_data.csv"

mv $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG/discordant_data.csv \
$CORE_PATH/$PROJECT_SAMPLE/REPORTS/CONCORDANCE_MS/$SM_TAG"_discordant_calls.csv"
