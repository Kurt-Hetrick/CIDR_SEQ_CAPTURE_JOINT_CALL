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
	HAPMAP_VCF=$4
	OMNI_VCF=$5
	ONEKG_SNPS_VCF=$6
	DBSNP_138_VCF=$7

	CORE_PATH=$8
	PROJECT_MS=$9
	PREFIX=${10}
	R_DIRECTORY=${11}
		export PATH=.:$R_DIRECTORY:$PATH
	SEND_TO=${12}

	# explicitly state the maximum number of gaussians to start of with

		MAX_GAUSSIANS="8"

START_VQSR_SNV=`date '+%s'`

	CMD=$JAVA_1_8'/java -jar'
	CMD=$CMD' '$GATK_DIR'/GenomeAnalysisTK.jar'
	CMD=$CMD' -T VariantRecalibrator'
	CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'
	CMD=$CMD' -R '$REF_GENOME
	CMD=$CMD' --input:VCF '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.raw.vcf'
	CMD=$CMD' -resource:hapmap,known=false,training=true,truth=true,prior=15.0 '$HAPMAP_VCF
	CMD=$CMD' -resource:omni,known=false,training=true,truth=true,prior=12.0 '$OMNI_VCF
	CMD=$CMD' -resource:1000G,known=false,training=true,truth=false,prior=10.0 '$ONEKG_SNPS_VCF
	CMD=$CMD' -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 '$DBSNP_138_VCF
	CMD=$CMD' -recalFile '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.SNP.recal'
	CMD=$CMD' -tranchesFile '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.SNP.tranches'
	CMD=$CMD' -rscriptFile '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.SNP.R'
	CMD=$CMD' -tranche 100.0'
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
	CMD=$CMD' -mode SNP'
	CMD=$CMD' -an QD'
	CMD=$CMD' -an MQRankSum'
	CMD=$CMD' -an MQ'
	CMD=$CMD' -an ReadPosRankSum'
	CMD=$CMD' -an SOR'
	CMD=$CMD' -an FS'
	CMD=$CMD' --maxGaussians '$MAX_GAUSSIANS

	echo $CMD | bash

	# capture the exit status

		SCRIPT_STATUS=`echo $?`

	# if vqsr fails then retry by decrementing the number of max gaussians by 1 until you get to 4
	# if it still does not work after setting it to four then stop trying

		if [ $SCRIPT_STATUS -ne 0 ]
			then
				until [[ $SCRIPT_STATUS -eq 0 || $MAX_GAUSSIANS -le 4 ]]
					do
						CMD=$(echo $CMD | sed 's/ --maxGaussians '"$MAX_GAUSSIANS"'//g')
						MAX_GAUSSIANS=$[$MAX_GAUSSIANS-1]
						CMD=$CMD' --maxGaussians '$MAX_GAUSSIANS
						echo $CMD | bash
						SCRIPT_STATUS=`echo $?`
				done
		fi

	# if it fails the first time but ultimately works send a notification saying that the parameter has changed and that the methods document needs to change for release
	# if it ends up failing altogether send a notification saying that I need to look at it.
	# will probably have to start with removing the MQ annotation and go from there.

		if [[ $SCRIPT_STATUS -eq 0 && $MAX_GAUSSIANS -ge 4 && $MAX_GAUSSIANS -lt 8 ]]
			then
				printf "The number of max Gaussians has been changed to $MAX_GAUSSIANS instead of 8 for\n \
				PROJECT:\n \
				$PROJECT_MS\n \
				VCF PREFIX:\n \
				$PREFIX\n \
				This needs to be reflected in the methods release document." \
				| mail -s "SNP VariantRecalibrator parameter changed for $PROJECT_MS" \
				$SEND_TO
			elif [ $SCRIPT_STATUS -ne 0 ]
				then
					printf "This has failed SNP VariantRecalibrator and Kurt needs to look at this for:\n \
					PROJECT:\n \
					$PROJECT_MS\n \
					VCF PREFIX:\n \
					$PREFIX" \
					| mail -s "SNP VariantRecalibrator FAILED for $PROJECT_MS" \
					$SEND_TO
			else
			:
		fi

echo $CMD >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"

END_VQSR_SNV=`date '+%s'`

echo $PROJECT_MS",E01,VQSR_SNV,"$HOSTNAME","$START_VQSR_SNV","$END_VQSR_SNV \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"

# write out SNP MAX GAUSSIAN SETTING TO FILE

	echo "SNP MAX GAUSSIANS SETTING FOR ${PROJECT_MS} at $(date) :::: ${MAX_GAUSSIANS}" \
	>> ${CORE_PATH}/${PROJECT_MS}/${PROJECT_MS}.VQSR_SETTINGS.txt

# exit with the signal from the program

	exit $SCRIPT_STATUS
