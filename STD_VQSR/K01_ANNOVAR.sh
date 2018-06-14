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

# INPUT PARAMETERS

CIDRSEQSUITE_ANNOVAR_JAVA=$1
CIDRSEQSUITE_DIR_4_0=$2
CIDRSEQSUITE_PROPS=$3

CORE_PATH=$4
PROJECT_MS=$5
PREFIX=$6

# I'm not sure why he is switching b/w different folders for all of these intermediate steps
# I might want to revisit and clean this up later

START_ANNOVAR=`date '+%s'`

zcat $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf.gz" \
$CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX".BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf"

# fyi after -annovar_directory_annotation the arguments are for
# 1. path to vcf file directory
# 2. path to output directory

$CIDRSEQSUITE_ANNOVAR_JAVA/java -jar \
-Duser.home=$CIDRSEQSUITE_PROPS \
-Xmx300g \
$CIDRSEQSUITE_DIR_4_0/CIDRSeqSuite.jar \
-pipeline \
-annovar_directory_annotation \
$CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX \
$CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX

echo \
$CIDRSEQSUITE_ANNOVAR_JAVA/java -jar \
-Duser.home=$CIDRSEQSUITE_PROPS \
-Xmx300g \
$CIDRSEQSUITE_DIR_4_0/CIDRSeqSuite.jar \
-pipeline \
-annovar_directory_annotation \
$CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX \
$CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX \
>> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"

echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"

du -ah $CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX \
	| egrep "csv|txt" \
	| awk '{print "mv -v",$2,"'$CORE_PATH'""/""'$PROJECT_MS'""/REPORTS/ANNOVAR/"}' \
	| bash

rm -rvf $CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX

END_ANNOVAR=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_MS",K01,ANNOVAR,"$HOSTNAME","$START_ANNOVAR","$END_ANNOVAR \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"

# check to see if the index is generated which should send an non-zero exit signal if not.
# eventually, will want to check the exit signal above and push out whatever it is at the end. Not doing that today though.

ls $CORE_PATH/$PROJECT_MS/REPORTS/ANNOVAR/$PREFIX".BEDsuperset.VQSR.1KG.ExAC3.REFINED.txt"
