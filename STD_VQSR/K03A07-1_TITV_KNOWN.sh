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

SAMTOOLS_0118_DIR=$1

CORE_PATH=$2
PROJECT_SAMPLE=$3
SM_TAG=$4

START_KNOWN_TITV=`date '+%s'`

CMD=$SAMTOOLS_0118_DIR'/bcftools/vcfutils.pl qstats '$CORE_PATH'/'$PROJECT_SAMPLE'/TEMP/'$SM_TAG'.Release.Known.TiTv.vcf'
CMD=$CMD' >| '$CORE_PATH'/'$PROJECT_SAMPLE'/REPORTS/TI_TV_MS/'$SM_TAG'_Known_.titv.txt'

echo $CMD >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
echo $CMD | bash

END_KNOWN_TITV=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_SAMPLE",M01,KNOWN_TITV,"$HOSTNAME","$START_KNOWN_TITV","$END_KNOWN_TITV \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".JOINT.CALL.WALL.CLOCK.TIMES.csv"
