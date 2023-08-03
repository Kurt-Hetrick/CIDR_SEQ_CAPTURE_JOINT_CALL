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

	echo

# INPUT PARAMETERS

	JAVA_1_8=$1
	GATK_DIR=$2
	REF_GENOME=$3
	CORE_PATH=$4

	PROJECT_MS=$5
	PREFIX=$6
	SNP_SENSITIVITY=$7

START_APPLY_VQSR_SNV=$(date '+%s') # capture time process starts for wall clock tracking purposes.

# construct command line

	CMD="${JAVA_1_8}/java -jar"
	CMD=${CMD}" ${GATK_DIR}/GenomeAnalysisTK.jar"
	CMD=${CMD}" -T ApplyRecalibration"
	CMD=${CMD}" --disable_auto_index_creation_and_locking_when_reading_rods"
	CMD=${CMD}" -R ${REF_GENOME}"
	CMD=${CMD}" --input:VCF ${CORE_PATH}/${PROJECT_MS}/MULTI_SAMPLE/${PREFIX}.raw.vcf"
	CMD=${CMD}" -recalFile ${CORE_PATH}/${PROJECT_MS}/MULTI_SAMPLE/${PREFIX}.SNP.recal"
	CMD=${CMD}" -tranchesFile ${CORE_PATH}/${PROJECT_MS}/MULTI_SAMPLE/${PREFIX}.SNP.tranches"
	CMD=${CMD}" -o ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.SNP.VQSR.vcf"
	CMD=${CMD}" -mode SNP"
	CMD=${CMD}" --ts_filter_level ${SNP_SENSITIVITY}"

	# write command line to file and execute the command line

		echo ${CMD} >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${PROJECT_MS}_command_lines.txt
		echo >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${PROJECT_MS}_command_lines.txt
		echo ${CMD} | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=$(echo $?)

		# if exit does not equal 0 then exit with whatever the exit signal is at the end.
		# also write to file that this job failed

			# if
			# 	[ "${SCRIPT_STATUS}" -ne 0 ]
			# then
			# 	echo ${PROJECT_MS} ${HOSTNAME} ${JOB_NAME} ${USER} ${SCRIPT_STATUS} ${SGE_STDERR_PATH} \
			# 	>> ${CORE_PATH}/${PROJECT_MS}/TEMP/${SAMPLE_SHEET_NAME}_${SUBMIT_STAMP}_ERRORS.txt

			# 	exit ${SCRIPT_STATUS}
			# fi

END_APPLY_VQSR_SNV=$(date '+%s') # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${PROJECT_MS},F01,APPLY_VQSR_SNV,${HOSTNAME},${START_APPLY_VQSR_SNV},${END_APPLY_VQSR_SNV} \
	>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}.JOINT.CALL.WALL.CLOCK.TIMES.csv

# write out VQSR TRUTH SENSTIVITY SETTING TO FILE

	echo "VQSR SNP TRUTH SENSITIVITY SETTING FOR ${PROJECT_MS} AT $(date) :::: ${SNP_SENSITIVITY}%" \
	>> ${CORE_PATH}/${PROJECT_MS}/${PROJECT_MS}.VQSR_SETTINGS.txt

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
