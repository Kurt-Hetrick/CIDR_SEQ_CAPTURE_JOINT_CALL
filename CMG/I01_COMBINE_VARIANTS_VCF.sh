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
	GATK_DIR=$2
	REF_GENOME=$3
	CORE_PATH=$4
	
	PROJECT_MS=$5
	PREFIX=$6

START_COMBINE_VARIANTS=`date '+%s'`

	CMD=$JAVA_1_8'/java'
	CMD=$CMD' -jar '$GATK_DIR'/GenomeAnalysisTK.jar'
	CMD=$CMD' -T CombineVariants'
	CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'	
	CMD=$CMD' -R '$REF_GENOME
	CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT_MS'/TEMP/'$PREFIX'.HC.SNP.INDEL.VQSR.COMMON.BIALLELIC.vcf.gz'
	CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT_MS'/TEMP/'$PREFIX'.HC.SNP.INDEL.VQSR.MULTIALLELIC.vcf.gz'

	for VCF in $(ls $CORE_PATH/$PROJECT_MS/TEMP/$PREFIX".HC.SNP.INDEL.VQSR.RARE.BIALLELIC.ANNOTATED."*.vcf.gz)
	do
	  CMD=$CMD' --variant '$VCF
	done

	CMD=$CMD' --genotypemergeoption UNSORTED'
	CMD=$CMD' -o '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.vcf.gz'

	echo $CMD >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
	echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
	echo $CMD | bash

END_COMBINE_VARIANTS=`date '+%s'`

echo $PROJECT_MS",J01,COMBINE_VARIANTS,"$HOSTNAME","$START_COMBINE_VARIANTS","$END_COMBINE_VARIANTS \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"

# check to see if the index is generated which should send an non-zero exit signal if not.
# eventually, will want to check the exit signal above and push out whatever it is at the end. Not doing that today though.

ls $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.vcf.gz.tbi"
