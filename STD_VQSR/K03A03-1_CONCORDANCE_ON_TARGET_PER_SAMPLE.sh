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
# redirecting stderr/stdout to file as a log.

	set

	echo

# INPUT VARIABLES

	BEDTOOLS_DIR=$1
	JAVA_1_8=$2
	GATK_DIR_4011=$3
	CIDRSEQSUITE_7_5_0_DIR=$4

	CORE_PATH=$5
	PROJECT_SAMPLE=$6
	SM_TAG=$7
	PROJECT_MS=$8
	TARGET_BED=$9
		TARGET_BED_NAME=$(basename ${TARGET_BED} .bed)
	REF_GENOME=${10}
	PREFIX=${11}
	VERACODE_CSV=${12}

START_CONCORDANCE=`date '+%s'` # capture time process starts for wall clock tracking purposes.

# look for a final report and store it as a variable. if there are multiple ones, then take the newest one

	FINAL_REPORT_FILE_TEST=$(ls -tr ${CORE_PATH}/${PROJECT_SAMPLE}/Pretesting/Final_Genotyping_Reports/*${SM_TAG}* \
		| tail -n 1)

# if final report exists containing the full sm-tag, then cidrseqsuite magic

	if [[ ! -z "${FINAL_REPORT_FILE_TEST}" ]]
		then
			FINAL_REPORT=${FINAL_REPORT_FILE_TEST}

	# if it does not exist, and if the $SM_TAG does not begin with an integer then split $SM_TAG On a @ or _ or -
	# look for a final report that contains that that first element of the $SM_TAG
	# bonus feature. if this first tests true but the file still does not exist then cidrseqsuite magic files b/c no file exists

	elif [[ ${SM_TAG} != [0-9]* ]]
		then
			# note that underscore has to be before dash in bash regular expression
			HAPMAP=${SM_TAG%%[@_-]*}

			FINAL_REPORT=$(ls ${CORE_PATH}/${PROJECT_SAMPLE}/Pretesting/Final_Genotyping_Reports/*${HAPMAP}* \
				| head -n 1)

			# if there is no report for a hapmap sample then exit program with code 1

			if [[ -z "${FINAL_REPORT}" ]]
			then

				echo
				echo At this time, you are looking for a final report that does not exist or fails to meet the current logic for finding a final report.
				echo Please talk to Kurt, because he loves to talk.
				echo

				FINAL_REPORT="FILE_DOES_NOT_EXIST"
				exit 1

			fi

	else

	# both conditions fails then echo the below message and give a dummy value for the $FINAL_REPORT
	# exit with 1

		echo
		echo At this time, you are looking for a final report that does not exist or fails to meet the current logic for finding a final report.
		echo Please talk to Kurt, because he loves to talk.
		echo

		FINAL_REPORT="FILE_DOES_NOT_EXIST"
		exit 1

	fi

# find what row the metadata header ends on.
# then add 1 to find where where the column/field headers start on.

	METADATA_HEADER_END_ROW=$(egrep -m 1 -n "^\[Data\]" ${FINAL_REPORT} \
		| cut -d ":" -f 1)

	FIELD_HEADER_ROW=$(echo ${METADATA_HEADER_END_ROW} \
		| awk '{print $0+1}')

# find what field number the CHR, POS, SNP NAME, SNP field columns are on.
# CHR AND POS (1-based) are used to generate bed file.
# SNP NAME isn't technically necessary, but could be nice to have for investigative purposes.
# SNP field is to remove insertion/deletions...can't remember what the copy number designation is.

	CHR_FIELD_NUMBER=(`sed 's/\r//g' ${FINAL_REPORT} \
		| awk 'NR=="'$FIELD_HEADER_ROW'"' \
		| sed 's/,/\n/g' \
		| cat -n \
		| sed 's/^ *//g' \
		| awk 'BEGIN {FS="\t"} \
			$2=="Chr" \
			{print $1}'`)

	POSITION_FIELD_NUMBER=(`sed 's/\r//g' ${FINAL_REPORT} \
		| awk 'NR=="'$FIELD_HEADER_ROW'"' \
		| sed 's/,/\n/g' \
		| cat -n \
		| sed 's/^ *//g' \
		| awk 'BEGIN {FS="\t"} \
			$2=="Position" \
			{print $1}'`)

	SNP_NAME_FIELD_NUMBER=(`sed 's/\r//g' ${FINAL_REPORT} \
		| awk 'NR=="'$FIELD_HEADER_ROW'"' \
		| sed 's/,/\n/g' \
		| cat -n \
		| sed 's/^ *//g' \
		| awk 'BEGIN {FS="\t"} \
			$2=="SNP Name" \
			{print $1}'`)

	SNP_FIELD_NUMBER=(`sed 's/\r//g' ${FINAL_REPORT} \
		| awk 'NR=="'$FIELD_HEADER_ROW'"' \
		| sed 's/,/\n/g' \
		| cat -n \
		| sed 's/^ *//g' \
		| awk 'BEGIN {FS="\t"} \
			$2=="SNP" \
			{print $1}'`)

# make a bed file from final report removing indel and only using karyotype chromosomes
# intersect with the (fixed) target bed file.

	sed 's/\r//g' ${FINAL_REPORT} \
		| awk -v FIELD_HEADER_ROW="$FIELD_HEADER_ROW" \
			'NR>FIELD_HEADER_ROW' \
		| awk \
			-v CHR_FIELD_NUMBER="$CHR_FIELD_NUMBER" \
			-v POSITION_FIELD_NUMBER="$POSITION_FIELD_NUMBER" \
			-v SNP_NAME_FIELD_NUMBER="$SNP_NAME_FIELD_NUMBER" \
			-v SNP_FIELD_NUMBER="$SNP_FIELD_NUMBER" \
			'BEGIN {FS=",";OFS="\t"} \
					$SNP_FIELD_NUMBER!~"I" \
					{print $CHR_FIELD_NUMBER,\
					$POSITION_FIELD_NUMBER-1,\
					$POSITION_FIELD_NUMBER,\
					$SNP_NAME_FIELD_NUMBER}' \
		| egrep "^[1-9]|^X|^Y" \
		| $BEDTOOLS_DIR/bedtools intersect \
			-a - \
			-b ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}.bed \
	>| ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}_FINAL_REPORT.bed

# Extract out sample, remove non-passing, non-variant from final ms vcf.
# using only those positions in the on target final report bed file.

	CMD="${JAVA_1_8}/java -jar"
		CMD=${CMD}" ${GATK_DIR_4011}/gatk-package-4.0.11.0-local.jar"
	CMD=${CMD}" SelectVariants"
		CMD=${CMD}" --reference ${REF_GENOME}"
		CMD=${CMD}" --variant ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf"
		CMD=${CMD}" --exclude-non-variants"
		CMD=${CMD}" --exclude-filtered"
		CMD=${CMD}" --remove-unused-alternates"
		CMD=${CMD}" --keep-original-ac"
		CMD=${CMD}" -L ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}_FINAL_REPORT.bed"
		CMD=${CMD}" --sample-name ${SM_TAG}"
	CMD=${CMD}" --output ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}_MS_FINAL_REPORT_ON_TARGET.vcf"

# write command line to file and execute the command line

	echo ${CMD} >> ${CORE_PATH}/${PROJECT_SAMPLE}/COMMAND_LINES/${SM_TAG}_command_lines.txt
	echo >> ${CORE_PATH}/${PROJECT_SAMPLE}/COMMAND_LINES/${SM_TAG}_command_lines.txt
	echo ${CMD} | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

		# if exit does not equal 0 then exit with whatever the exit signal is at the end.
		# also write to file that this job failed

			if [ "${SCRIPT_STATUS}" -ne 0 ]
				then
					echo ${SM_TAG} ${HOSTNAME} ${JOB_NAME} ${USER} ${SCRIPT_STATUS} ${SGE_STDERR_PATH} \
					>> ${CORE_PATH}/${PROJECT_MS}/TEMP/${SAMPLE_SHEET_NAME}_${SUBMIT_STAMP}_ERRORS.txt
					exit ${SCRIPT_STATUS}
			fi

# -single_sample_concordance
# Performs concordance between one vcf file and one final report. The vcf must be single sample.

	CMD="${JAVA_1_8}/java -jar"
		CMD=${CMD}" ${CIDRSEQSUITE_7_5_0_DIR}/CIDRSeqSuite.jar"
		CMD=${CMD}" -single_sample_concordance"
		# [1] path_to_vcf_file
		CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}_MS_FINAL_REPORT_ON_TARGET.vcf"
		# [2] path_to_final_report_file
		CMD=${CMD}" ${FINAL_REPORT}"
		# [3] path_to_bed_file
		CMD=${CMD}" ${TARGET_BED}"
		# [4] path_to_liftover_file
		CMD=${CMD}" ${VERACODE_CSV}"
	# [5] path_to_output_directory
	CMD=${CMD}" ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/CONCORDANCE_MS/"

# write command line to file and execute the command line

	echo ${CMD} >> ${CORE_PATH}/${PROJECT_SAMPLE}/COMMAND_LINES/${SM_TAG}_command_lines.txt
	echo >> ${CORE_PATH}/${PROJECT_SAMPLE}/COMMAND_LINES/${SM_TAG}_command_lines.txt
	echo ${CMD} | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

		# if exit does not equal 0 then exit with whatever the exit signal is at the end.
		# also write to file that this job failed

			if [ "${SCRIPT_STATUS}" -ne 0 ]
				then
					echo ${SM_TAG} ${HOSTNAME} ${JOB_NAME} ${USER} ${SCRIPT_STATUS} ${SGE_STDERR_PATH} \
					>> ${CORE_PATH}/${PROJECT_MS}/TEMP/${SAMPLE_SHEET_NAME}_${SUBMIT_STAMP}_ERRORS.txt
					exit ${SCRIPT_STATUS}
			fi

END_CONCORDANCE=`date '+%s'` # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${PROJECT_SAMPLE},M01,CONCORDANCE,${HOSTNAME},${START_CONCORDANCE},${END_CONCORDANCE} \
	>> ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/${PROJECT_SAMPLE}.JOINT.CALL.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
