# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge that you are in the users current working directory
#$ -cwd

# tell sge to export the users environment variables
#$ -V

# tell sge to submit at this priority setting
#$ -p -10

# tell sge to output both stderr and stdout to the same file
#$ -j y

# export all variables, useful to find out what compute node the program was executed on

	set

	echo

# INPUT VARIABLES

	JAVA_1_8=$1
	GATK_DIR_4011=$2

	REF_DICT=$3
	KNOWN_SNPS=$4
	CORE_PATH=$5
	PROJECT_MS=$6
	PROJECT_TITV_BED=$7
		PROJECT_TITV_BED_NAME=$(basename ${PROJECT_TITV_BED} .bed)
	PREFIX=$8
	SAMPLE_SHEET=$9
		SAMPLE_SHEET_NAME=$(basename ${SAMPLE_SHEET} .csv)
	SUBMIT_STAMP=${10}

# FILTER INDEL AND MIXED VARIANTS

START_VCF_METRICS_TITV=$(date '+%s') # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="${JAVA_1_8}/java -jar"
			CMD=${CMD}" ${GATK_DIR_4011}/gatk-package-4.0.11.0-local.jar"
		CMD=${CMD}" CollectVariantCallingMetrics"
			CMD=${CMD}" --SEQUENCE_DICTIONARY ${REF_DICT}"
			CMD=${CMD}" --DBSNP ${KNOWN_SNPS}"
			CMD=${CMD}" --THREAD_COUNT 4"
			CMD=${CMD}" --TARGET_INTERVALS ${CORE_PATH}/${PROJECT_MS}/TEMP/${PROJECT_MS}-${PROJECT_TITV_BED_NAME}-picard.bed"
			CMD=${CMD}" --INPUT ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.FILTERED.GT.REFINED.vcf"
		CMD=${CMD}" --OUTPUT ${CORE_PATH}/${PROJECT_MS}/REPORTS/VCF_METRICS/MULTI_SAMPLE/${PREFIX}_TITV"
		CMD=${CMD}" &&"
		CMD=${CMD}" mv -v"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/REPORTS/VCF_METRICS/MULTI_SAMPLE/${PREFIX}_TITV.variant_calling_detail_metrics"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/REPORTS/VCF_METRICS/MULTI_SAMPLE/${PREFIX}_TITV.variant_calling_detail_metrics.txt"
		CMD=${CMD}" &&"
		CMD=${CMD}" mv -v"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/REPORTS/VCF_METRICS/MULTI_SAMPLE/${PREFIX}_TITV.variant_calling_summary_metrics"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/REPORTS/VCF_METRICS/MULTI_SAMPLE/${PREFIX}_TITV.variant_calling_summary_metrics.txt"

	# write command line to file and execute the command line

		echo ${CMD} >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${PROJECT_MS}_command_lines.txt
		echo >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${PROJECT_MS}_command_lines.txt
		echo ${CMD} | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=$(echo $?)

		# if exit does not equal 0 then exit with whatever the exit signal is at the end.
		# also write to file that this job failed

			if
				[ "${SCRIPT_STATUS}" -ne 0 ]
			then
				echo ${PROJECT_MS} ${HOSTNAME} ${JOB_NAME} ${USER} ${SCRIPT_STATUS} ${SGE_STDERR_PATH} \
				>> ${CORE_PATH}/${PROJECT_MS}/TEMP/${SAMPLE_SHEET_NAME}_${SUBMIT_STAMP}_ERRORS.txt
				exit ${SCRIPT_STATUS}
			fi

END_VCF_METRICS_TITV=$(date '+%s') # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${PROJECT_MS},K01,VCF_METRICS_TITV,${HOSTNAME},${START_VCF_METRICS_TITV},${END_VCF_METRICS_TITV} \
	>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
