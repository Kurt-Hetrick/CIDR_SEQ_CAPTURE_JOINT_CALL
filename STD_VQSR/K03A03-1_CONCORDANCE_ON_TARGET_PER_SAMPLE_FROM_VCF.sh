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

	GENOTYPING_VCF_TEST=$(ls -tr ${CORE_PATH}/${PROJECT_SAMPLE}/Pretesting/Final_Genotyping_VCFs/*${SM_TAG}* \
		| tail -n 1)

# if final report exists containing the full sm-tag, then cidrseqsuite magic

	if [[ ! -z "${GENOTYPING_VCF_TEST}" ]]
		then
			GENOTYPING_VCF=${GENOTYPING_VCF_TEST}

	# if it does not exist, and if the $SM_TAG does not begin with an integer then split $SM_TAG On a @ or _ or -
	# look for a final report that contains that that first element of the $SM_TAG
	# bonus feature. if this first tests true but the file still does not exist then cidrseqsuite magic files b/c no file exists

	elif [[ ${SM_TAG} != [0-9]* ]]
		then
			# note that underscore has to be before dash in bash regular expression
			HAPMAP=${SM_TAG%%[@_-]*}

			GENOTYPING_VCF=$(ls ${CORE_PATH}/${PROJECT_SAMPLE}/Pretesting/Final_Genotyping_VCFs/*${HAPMAP}* \
				| head -n 1)

			# if there is no report for a hapmap sample then exit program with code 1

			if [[ -z "${GENOTYPING_VCF}" ]]
			then

				echo
				echo At this time, you are looking for a final report that does not exist or fails to meet the current logic for finding a final report.
				echo Please talk to Kurt, because he loves to talk.
				echo

				GENOTYPING_VCF="FILE_DOES_NOT_EXIST"
				exit 1

			fi

	else

	# both conditions fails then echo the below message and give a dummy value for the $GENOTYPING_VCF
	# exit with 1

		echo
		echo At this time, you are looking for a final report that does not exist or fails to meet the current logic for finding a final report.
		echo Please talk to Kurt, because he loves to talk.
		echo

		GENOTYPING_VCF="FILE_DOES_NOT_EXIST"
		exit 1

	fi

# intersect GENOTYPING VCF with the (fixed) target bed file.

	${BEDTOOLS_DIR}/bedtools intersect \
		-header \
		-a ${GENOTYPING_VCF} \
		-b ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}.bed \
	>| ${CORE_PATH}/${PROJECT_SAMPLE}/TEMP/${SM_TAG}/${SM_TAG}.TRUTH.TARGET_ONLY.vcf

# CONVERT ON TARGET GENOTYPING VCF INTO A BED FILE

	/mnt/linuxtools/BEDOPS/bin/vcf2bed \
		< ${CORE_PATH}/${PROJECT_SAMPLE}/TEMP/${SM_TAG}/${SM_TAG}.TRUTH.TARGET_ONLY.vcf \
	>| ${CORE_PATH}/${PROJECT_SAMPLE}/TEMP/${SM_TAG}/${SM_TAG}.TRUTH.TARGET_ONLY.bed

# CONVERT ON TARGET GENOTYPING BED FILE INTO PICARD INTERVAL FILE

	(grep ^@ ${REF_DICT} ; \
	cut -f 1-4 ${CORE_PATH}/${PROJECT_SAMPLE}/TEMP/${SM_TAG}/${SM_TAG}.TRUTH.TARGET_ONLY.bed \
		| awk 'BEGIN {OFS="\t"} \
			{print $1,$2,$3,"+",$4}' ) \
	>| ${CORE_PATH}/${PROJECT_SAMPLE}/TEMP/${SM_TAG}/${SM_TAG}.TRUTH.TARGET_ONLY.picard.bed

# THE GTC TO VCF CONVERTOR USES THE THE BEADCHIP POSITION TO POPULATE THE SAMPLE NAME IN THE VCF
# NEED TO OBTAIN TO MAP THIS NAME TO THE NAME IN THE MULTI-SAMPLE VCF
# I THINK SINCE THEY DON'T MATCH YOU HAVE TO SPECIFY THAT THEY ARE SUPPOSED TO MATCH

	BEADCHIP_POSITION=$(zgrep -m 1 "^#CHROM" $GENOTYPING_VCF | cut -f 10)

# RUN GENOTYPING CONCORDANCE

	CMD="${JAVA_1_8}/java -jar"
		CMD=${CMD}" /mnt/linuxtools/PICARD/picard-2.26.2/picard.jar"
	CMD=${CMD}" GenotypeConcordance"
		CMD=${CMD}" CALL_VCF=${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.FILTERED.GT.REFINED.vcf"
		CMD=${CMD}" CALL_SAMPLE=${SM_TAG}"
		CMD=${CMD}" TRUTH_VCF=${GENOTYPING_VCF}"
		CMD=${CMD}" TRUTH_SAMPLE=${BEADCHIP_POSITION}"
		CMD=${CMD}" INTERVALS=${CORE_PATH}/${PROJECT_SAMPLE}/TEMP/${SM_TAG}/${SM_TAG}.TRUTH.TARGET_ONLY.picard.bed"
		CMD=${CMD}" USE_VCF_INDEX=TRUE"
		CMD=${CMD}" OUTPUT_VCF=TRUE"
	CMD=${CMD}" OUTPUT=${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/CONCORDANCE_METRICS/RELEASE/${SM_TAG}"

# write command line to file and execute the command line

	echo ${CMD} >> ${CORE_PATH}/${PROJECT_SAMPLE}/COMMAND_LINES/${SM_TAG}.COMMAND.LINES.txt
	echo >> ${CORE_PATH}/${PROJECT_SAMPLE}/COMMAND_LINES/${SM_TAG}.COMMAND.LINES.txt
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

# write command line to file and execute the command line

	echo ${CMD} >> ${CORE_PATH}/${PROJECT_SAMPLE}/COMMAND_LINES/${SM_TAG}.COMMAND.LINES.txt
	echo >> ${CORE_PATH}/${PROJECT_SAMPLE}/COMMAND_LINES/${SM_TAG}.COMMAND.LINES.txt
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
