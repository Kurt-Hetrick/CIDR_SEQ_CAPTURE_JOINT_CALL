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
	PROJECT_MS=$4
	PREFIX=$5
	GRCH38_REF=$6
	HG19_TO_GRCH38_CHAIN=$7

# liftover from hg19 to GRCh38

START_LIFTOVER_MS_HG38=`date '+%s'`

	$JAVA_1_8/java -jar \
	$PICARD_DIR/picard.jar \
	LiftoverVcf \
	INPUT=$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.HG19.LIFTOVER.vcf.gz" \
	OUTPUT=$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.GRCh38.LIFTOVER.vcf.gz" \
	REJECT=$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.GRCh38.LIFTOVER.REJECTED.vcf.gz" \
	REFERENCE_SEQUENCE=$GRCH38_REF \
	CHAIN=$HG19_TO_GRCH38_CHAIN

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

END_LIFTOVER_MS_HG38=`date '+%s'`

echo $SM_TAG"_"$PROJECT_MS",L.01,LIFTOVER_REFINED_VCF_HG38,"$HOSTNAME","$START_LIFTOVER_MS_HG38","$END_LIFTOVER_MS_HG38 \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".WALL.CLOCK.TIMES.csv"

echo \
$JAVA_1_8/java -jar \
$PICARD_DIR/picard.jar \
LiftoverVcf \
INPUT=$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.HG19.LIFTOVER.vcf.gz" \
OUTPUT=$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.GRCh38.LIFTOVER.vcf.gz" \
REJECT=$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.GRCh38.LIFTOVER.REJECTED.vcf.gz" \
REFERENCE_SEQUENCE=$GRCH38_REF \
CHAIN=$HG19_TO_GRCH38_CHAIN \
>> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"-"$PREFIX".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"-"$PREFIX".COMMAND.LINES.txt"

# exit with the signal from the program

	exit $SCRIPT_STATUS
