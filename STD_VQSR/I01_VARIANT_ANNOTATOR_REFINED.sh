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
	DBSNP=$4

	CORE_PATH=$5
	PROJECT_MS=$6
	PREFIX=$7
	BED_FILE_NAME=$8

START_VARIANT_ANNOTATOR_2=$(date '+%s') # capture time process starts for wall clock tracking purposes.

# construct cmd line

	CMD="${JAVA_1_8}/java -jar"
	CMD=${CMD}" ${GATK_DIR}/GenomeAnalysisTK.jar"
	CMD=${CMD}" -T VariantAnnotator"
	CMD=${CMD}" --disable_auto_index_creation_and_locking_when_reading_rods"
	CMD=${CMD}" -R ${REF_GENOME}"
	CMD=${CMD}" --variant ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.${BED_FILE_NAME}.GT.REFINED.temp.vcf.gz"
	CMD=${CMD}" -o ${CORE_PATH}/${PROJECT_MS}/TEMP/${BED_FILE_NAME}.r.vcf.gz"
	CMD=${CMD}" -L ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.${BED_FILE_NAME}.GT.REFINED.temp.vcf.gz"
	CMD=${CMD}" -A GenotypeSummaries"

# write command line to file and execute the command line

	echo ${CMD} >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${PROJECT_MS}_command_lines.txt
	echo >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${PROJECT_MS}_command_lines.txt
	echo ${CMD} | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=$(echo $?)

		### currently not being implemented.
		# if exit does not equal 0 then exit with whatever the exit signal is at the end.
		# also write to file that this job failed.

			# if
			# 	[ "${SCRIPT_STATUS}" -ne 0 ]
			# then
			# 	echo ${SM_TAG} ${HOSTNAME} ${JOB_NAME} ${USER} ${SCRIPT_STATUS} ${SGE_STDERR_PATH} \
			# 	>> ${CORE_PATH}/${PROJECT}/TEMP/${SAMPLE_SHEET_NAME}_${SUBMIT_STAMP}_ERRORS.txt
			# 	exit ${SCRIPT_STATUS}
			# fi

END_VARIANT_ANNOTATOR_2=$(date '+%s') # capture time process stops for wall clock tracking purposes.

# write wall clock times to file

	echo ${PROJECT_MS},I01,VARIANT_ANNOTATOR_2,${HOSTNAME},${START_VARIANT_ANNOTATOR_2},${END_VARIANT_ANNOTATOR_2} \
	>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}.JOINT.CALL.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
