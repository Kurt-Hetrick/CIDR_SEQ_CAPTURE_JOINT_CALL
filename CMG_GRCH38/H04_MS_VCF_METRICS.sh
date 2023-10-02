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

# INPUT PARAMETERS

	JAVA_1_8=$1
	GATK_DIR_4011=$2
	CORE_PATH=$3

	PROJECT_MS=$4
	PREFIX=$5
	DBSNP_138_VCF=$6
	REF_DICT=$7

START_MS_VCF_METRICS=`date '+%s'`

	CMD=$JAVA_1_8'/java -jar'
	CMD=$CMD' '$GATK_DIR_4011'/gatk-package-4.0.11.0-local.jar'
	CMD=$CMD' CollectVariantCallingMetrics'
	CMD=$CMD' --INPUT '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.HC.SNP.INDEL.VQSR.vcf.gz'
	CMD=$CMD' --DBSNP '$DBSNP_138_VCF
	CMD=$CMD' --SEQUENCE_DICTIONARY '$REF_DICT
	CMD=$CMD' --OUTPUT '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX
	CMD=$CMD' --THREAD_COUNT 4'

START_MS_VCF_METRICS=`date '+%s'`

echo $CMD >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo $CMD | bash

# add text extension to files.

	mv -v $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".variant_calling_detail_metrics" \
	$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".variant_calling_detail_metrics.txt"

	mv -v $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".variant_calling_summary_metrics" \
	$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".variant_calling_summary_metrics.txt"

echo $PROJECT_MS",H01,MS_VCF_METRICS,"$HOSTNAME","$START_MS_VCF_METRICS","$END_MS_VCF_METRICS \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"
