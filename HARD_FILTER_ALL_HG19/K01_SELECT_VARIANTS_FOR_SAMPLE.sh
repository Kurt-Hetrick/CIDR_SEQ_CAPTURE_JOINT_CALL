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
	GATK_DIR_4011=$2
	REF_GENOME=$3

	CORE_PATH=$4
	PROJECT_SAMPLE=$5
	PROJECT_MS=$6
	SM_TAG=$7
	PREFIX=$8

START_SELECT_ALL_SAMPLE=`date '+%s'`

# Extract out sample, remove non-passing, non-variant

	CMD=$JAVA_1_8'/java -jar'
	CMD=$CMD' '$GATK_DIR_4011'/gatk-package-4.0.11.0-local.jar'
	CMD=$CMD' SelectVariants'
	CMD=$CMD' --reference '$REF_GENOME
	CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT_MS'/TEMP/'$PREFIX'.HF.1KG.ExAC3.REFINED.vcf'
	CMD=$CMD' --output '$CORE_PATH'/'$PROJECT_SAMPLE'/VCF/RELEASE/FILTERED_ON_BAIT/'$SM_TAG'_MS_OnBait.vcf.gz'
	CMD=$CMD' --exclude-non-variants'
	CMD=$CMD' --exclude-filtered'
	CMD=$CMD' --remove-unused-alternates'
	CMD=$CMD' --keep-original-ac'
	CMD=$CMD' --sample-name '$SM_TAG

		echo $CMD >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
		echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
		echo $CMD | bash

END_SELECT_ALL_SAMPLE=`date '+%s'`

echo $PROJECT_SAMPLE",K01,SELECT_ALL_SAMPLE,"$HOSTNAME","$START_SELECT_ALL_SAMPLE","$END_SELECT_ALL_SAMPLE \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".JOINT.CALL.WALL.CLOCK.TIMES.csv"
