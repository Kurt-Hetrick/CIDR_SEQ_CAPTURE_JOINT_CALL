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

START_LIFTOVER_PHENODB=`date '+%s'`

	$JAVA_1_8/java -jar \
	$PICARD_DIR/picard.jar \
	LiftoverVcf \
	INPUT=$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz" \
	OUTPUT=$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT_HG19/$SM_TAG"_MS_OnBait.hg19liftover.vcf.gz" \
	REJECT=$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT_HG19/$SM_TAG"_MS_OnBait.hg19.hg19liftover.reject.vcf.gz" \
	REFERENCE_SEQUENCE=$HG19_REF \
	CHAIN=$HG38_TO_HG19_CHAIN

START_LIFTOVER_PHENODB=`date '+%s'`

echo $SM_TAG"_"$PROJECT_SAMPLE",M.01,LIFTOVER_PHENODB,"$HOSTNAME","$START_LIFTOVER_PHENODB","$END_LIFTOVER_PHENODB \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".WALL.CLOCK.TIMES.csv"

echo \
$JAVA_1_8/java -jar \
$PICARD_DIR/picard.jar \
LiftoverVcf \
INPUT=$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz" \
OUTPUT=$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT_HG19/$SM_TAG"_MS_OnBait.hg19liftover.vcf.gz" \
REJECT=$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT_HG19/$SM_TAG"_MS_OnBait.hg19.hg19liftover.reject.vcf.gz" \
REFERENCE_SEQUENCE=$HG19_REF \
CHAIN=$HG38_TO_HG19_CHAIN \
>> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
