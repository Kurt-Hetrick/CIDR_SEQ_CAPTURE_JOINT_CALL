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

# create a blank lane b/w the output variables and the program logging output
echo

# INPUT VARIABLES

JAVA_1_8=$1
GATK_DIR=$2
REF_GENOME=$3
KNOWN_SNPS=$4

CORE_PATH=$5
PROJECT_SAMPLE=$6
SM_TAG=$7
TITV_BED=$8

START_KNOWN_SNV_TITV=`date '+%s'`

CMD=$JAVA_1_8'/java -jar'
CMD=$CMD' '$GATK_DIR'/GenomeAnalysisTK.jar'
CMD=$CMD' -T SelectVariants'
CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'
CMD=$CMD' -R '$REF_GENOME
CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT_SAMPLE'/VCF/RELEASE/FILTERED_ON_BAIT/'$SM_TAG'_MS_OnBait.vcf.gz'
CMD=$CMD' -o '$CORE_PATH'/'$PROJECT_SAMPLE'/TEMP/'$SM_TAG'.Release.Known.TiTv.vcf'
CMD=$CMD' -selectType SNP'
CMD=$CMD' --keepOriginalAC'
CMD=$CMD' --excludeFiltered'
CMD=$CMD' --excludeNonVariants'
CMD=$CMD' -L '$TITV_BED
CMD=$CMD' --concordance '$KNOWN_SNPS

echo $CMD >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
echo $CMD | bash

END_KNOWN_SNV_TITV=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_SAMPLE",L01,KNOWN_SNV_TITV,"$HOSTNAME","$START_KNOWN_SNV_TITV","$END_KNOWN_SNV_TITV \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".JOINT.CALL.WALL.CLOCK.TIMES.csv"
