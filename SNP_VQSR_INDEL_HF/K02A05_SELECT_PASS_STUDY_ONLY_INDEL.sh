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
	HAPMAP_SAMPLE_LIST=$7

START_STUDY_INDELS_PASS=`date '+%s'`

# for the study samples...via excluding the hapmap samples
# select passing INDEL,MIXED,MNP,SYMBOLIC sites that are only polymorphic in the hapmap samples
# REALLY WE ARE GOING FOR NON-SNP SITES HERE

	CMD=$JAVA_1_8'/java -jar'
	CMD=$CMD' '$GATK_DIR_4011'/gatk-package-4.0.11.0-local.jar'
	CMD=$CMD' SelectVariants'
	CMD=$CMD' --reference '$REF_GENOME
	CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT_MS'/TEMP/'$PREFIX'.BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.vcf'
	CMD=$CMD' --output '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/'$PREFIX'.BEDsuperset.VQSR.INDEL.STUDY.SAMPLES.PASS.vcf'
	CMD=$CMD' --select-type-to-include INDEL'
	CMD=$CMD' --select-type-to-include MIXED'
	CMD=$CMD' --select-type-to-include MNP'
	CMD=$CMD' --select-type-to-include SYMBOLIC'
	CMD=$CMD' --exclude-non-variants'
	CMD=$CMD' --exclude-filtered'
	CMD=$CMD' --remove-unused-alternates'
	CMD=$CMD' --exclude-sample-name '$HAPMAP_SAMPLE_LIST

echo $CMD >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo $CMD | bash

END_STUDY_INDELS_PASS=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_MS",J01,STUDY_INDELS_PASS,"$HOSTNAME","$START_STUDY_INDELS_PASS","$END_STUDY_INDELS_PASS \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"

# check to see if the index is generated which should send an non-zero exit signal if not.
# eventually, will want to check the exit signal above and push out whatever it is at the end. Not doing that today though.

ls $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/$PREFIX".BEDsuperset.VQSR.INDEL.STUDY.SAMPLES.PASS.vcf.idx"
