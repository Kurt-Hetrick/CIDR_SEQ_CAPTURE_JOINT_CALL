# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge that you are in the users current working directory
#$ -cwd

# tell sge to export the users environment variables
#$ -V

# tell sge to submit at this priority setting
#$ -p -10

# tell sge to output both stderr and stdout to the same file
#$ -j y

# export all variables, useful to find out what compute node the program was executed on

	set

	echo

# INPUT VARIABLES

	JAVA_1_8=$1
	PICARD_DIR=$2
	
	CORE_PATH=$3
	PROJECT_SAMPLE=$4
	SM_TAG=$5
	HG19_REF=$6
	HG38_TO_HG19_CHAIN=$7

# mkdir a directory in TEMP for the SM tag to decompress the target vcf file into

	# mkdir -p $CORE_PATH/$PROJECT/TEMP/$SM_TAG

# liftover from hg38 to hg19

START_LIFTOVER_SNV_TARGET_PASS=`date '+%s'`

	$JAVA_1_8/java -jar \
	$PICARD_DIR/picard.jar \
	LiftoverVcf \
	INPUT=$CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_SNV.vcf.gz" \
	OUTPUT=$CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET_HG19/$SM_TAG"_MS_OnTarget_SNV.hg19liftover.vcf.gz" \
	REJECT=$CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET_HG19/$SM_TAG"_MS_OnTarget_SNV.hg19.hg19liftover.reject.vcf.gz" \
	REFERENCE_SEQUENCE=$HG19_REF \
	CHAIN=$HG38_TO_HG19_CHAIN

END_LIFTOVER_SNV_TARGET_PASS=`date '+%s'`

echo $SM_TAG"_"$PROJECT_SAMPLE",M.01,SNV_TARGET_PASS_LIFTOVER,"$HOSTNAME","$START_SNV_TARGET_PASS_LIFTOVER","$END_SNV_TARGET_PASS_LIFTOVER \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".WALL.CLOCK.TIMES.csv"

echo \
$JAVA_1_8/java -jar \
$PICARD_DIR/picard.jar \
LiftoverVcf \
INPUT=$CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_SNV.vcf.gz" \
OUTPUT=$CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET_HG19/$SM_TAG/$SM_TAG"_MS_OnTarget_SNV.hg19liftover.vcf.gz" \
REJECT=$CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET_HG19/$SM_TAG"_MS_OnTarget_SNV.hg19.hg19liftover.reject.vcf.gz" \
REFERENCE_SEQUENCE=$HG19_REF \
CHAIN=$HG38_TO_HG19_CHAIN \
>> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
