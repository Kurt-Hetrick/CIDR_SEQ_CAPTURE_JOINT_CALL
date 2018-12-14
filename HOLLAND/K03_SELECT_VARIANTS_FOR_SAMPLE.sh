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

	CORE_PATH=$4
	PROJECT_SAMPLE=$5
	PROJECT_MS=$6
	SM_TAG=$7
	PREFIX=$8
	TABIX_DIR=$9

START_SELECT_ALL_SAMPLE=`date '+%s'`

# Extract out sample, remove non-passing, non-variant

	CMD=$JAVA_1_8'/java -jar'
	CMD=$CMD' '$GATK_DIR'/GenomeAnalysisTK.jar'
	CMD=$CMD' -T SelectVariants'
	CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'
	CMD=$CMD' -R '$REF_GENOME
	CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.vcf.gz'
	CMD=$CMD' -o '$CORE_PATH'/'$PROJECT_SAMPLE'/TEMP/'$SM_TAG'_MS_OnBait.TEMP.vcf.gz'
	CMD=$CMD' --keepOriginalAC'
	CMD=$CMD' -ef'
	CMD=$CMD' -env'
	CMD=$CMD' --removeUnusedAlternates'
	CMD=$CMD' -sn '$SM_TAG

echo $CMD >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
echo $CMD | bash

# Now remove the records where the alternate allele is just a *

	( zgrep "^#" $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_OnBait.TEMP.vcf.gz" ; \
		zgrep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_OnBait.TEMP.vcf.gz" | awk '$5!="*"' ) \
	$TABIX_DIR/bgzip -c /dev/stdin \
	>| $CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz"

# index the vcf file

	$TABIX_DIR/tabix -f -p vcf $CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz"

END_SELECT_ALL_SAMPLE=`date '+%s'`

echo $PROJECT_SAMPLE",K01,SELECT_ALL_SAMPLE,"$HOSTNAME","$START_SELECT_ALL_SAMPLE","$END_SELECT_ALL_SAMPLE \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".JOINT.CALL.WALL.CLOCK.TIMES.csv"
