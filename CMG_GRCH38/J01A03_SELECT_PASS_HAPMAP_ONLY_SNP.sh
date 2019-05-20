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
	PROJECT_MS=$5
	PREFIX=$6
	STUDY_SAMPLE_LIST=$7

START_HAPMAP_SNP_PASS=`date '+%s'`

# for the hapmap samples...via excluding the study samples
# select passing SNP sites that are only polymorphic in the hapmap samples

	CMD=$JAVA_1_8'/java -jar'
	CMD=$CMD' '$GATK_DIR_4011'/gatk-package-4.0.11.0-local.jar'
	CMD=$CMD' SelectVariants'
	CMD=$CMD' --reference '$REF_GENOME
	CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.vcf.gz'
	CMD=$CMD' --output '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/'$PREFIX'.BEDsuperset.VQSR.SNP.HAPMAP.SAMPLES.PASS.vcf'
	CMD=$CMD' --select-type-to-include SNP'
	CMD=$CMD' --exclude-non-variants'
	CMD=$CMD' --exclude-filtered'
	CMD=$CMD' --remove-unused-alternates'
	CMD=$CMD' --exclude-sample-name '$STUDY_SAMPLE_LIST

echo $CMD >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo $CMD | bash

END_HAPMAP_SNP_PASS=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_MS",J01,HAPMAP_SNP_PASS,"$HOSTNAME","$START_HAPMAP_SNP_PASS","$END_HAPMAP_SNP_PASS \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"

# check to see if the index is generated which should send an non-zero exit signal if not.
# eventually, will want to check the exit signal above and push out whatever it is at the end. Not doing that today though.

ls $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/$PREFIX".BEDsuperset.VQSR.SNP.HAPMAP.SAMPLES.PASS.vcf.idx"
