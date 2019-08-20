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
	ONEKG_INDELS_VCF=$4
	AXIOM_VCF=$5
	DBSNP_138_VCF=$6

	CORE_PATH=$7
	PROJECT_MS=$8
	PREFIX=$9
	R_DIRECTORY=${10}
		export PATH=.:$R_DIRECTORY:$PATH

START_VQSR_INDEL=`date '+%s'`

	CMD=$JAVA_1_8'/java -jar'
	CMD=$CMD' '$GATK_DIR'/GenomeAnalysisTK.jar'
	CMD=$CMD' -T VariantRecalibrator'
	CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'
	CMD=$CMD' -R '$REF_GENOME
	CMD=$CMD' --input:VCF '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.raw.HC.vcf'
	CMD=$CMD' -resource:mills,known=false,training=true,truth=true,prior=12.0 '$ONEKG_INDELS_VCF
	CMD=$CMD' -resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 '$AXIOM_VCF
	CMD=$CMD' -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 '$DBSNP_138_VCF
	CMD=$CMD' -recalFile '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.HC.INDEL.recal'
	CMD=$CMD' -tranchesFile '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.HC.INDEL.tranches'
	CMD=$CMD' -rscriptFile '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.HC.INDEL.R'
	CMD=$CMD' -tranche 100.0'
	CMD=$CMD' -tranche 99.99'
	CMD=$CMD' -tranche 99.98'
	CMD=$CMD' -tranche 99.97'
	CMD=$CMD' -tranche 99.96'
	CMD=$CMD' -tranche 99.95'
	CMD=$CMD' -tranche 99.94'
	CMD=$CMD' -tranche 99.93'
	CMD=$CMD' -tranche 99.92'
	CMD=$CMD' -tranche 99.91'
	CMD=$CMD' -tranche 99.9'
	CMD=$CMD' -tranche 99.8'
	CMD=$CMD' -tranche 99.7'
	CMD=$CMD' -tranche 99.6'
	CMD=$CMD' -tranche 99.5'
	CMD=$CMD' -tranche 99.4'
	CMD=$CMD' -tranche 99.3'
	CMD=$CMD' -tranche 99.2'
	CMD=$CMD' -tranche 99.1'
	CMD=$CMD' -tranche 99.0'
	CMD=$CMD' -tranche 98.0'
	CMD=$CMD' -tranche 97.0'
	CMD=$CMD' -tranche 96.0'
	CMD=$CMD' -tranche 95.0'
	CMD=$CMD' -tranche 90.0'
	CMD=$CMD' -mode INDEL'
	CMD=$CMD' --maxGaussians 4'
	CMD=$CMD' -an MQRankSum'
	CMD=$CMD' -an SOR'
	CMD=$CMD' -an ReadPosRankSum'
	CMD=$CMD' -an QD'
	CMD=$CMD' -an FS'

echo $CMD >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo $CMD | bash

END_VQSR_INDEL=`date '+%s'`

echo $PROJECT_MS",E01,VQSR_INDEL,"$HOSTNAME","$START_VQSR_INDEL","$END_VQSR_INDEL \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"

# check to see if the index is generated which should send an non-zero exit signal if not.
# eventually, will want to check the exit signal above and push out whatever it is at the end. Not doing that today though.
# this is placeholder here...have to see what is generated. I think I want to check if a *recal.idx is generated or not.
