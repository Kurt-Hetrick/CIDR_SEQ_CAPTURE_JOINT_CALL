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
	TABIX_DIR=$2
	GATK_DIR=$3
	GRCH38_REF=$4
	CORE_PATH=$5
	PROJECT_MS=$6
	PREFIX=$7

START_CAT_VARIANTS=$(date '+%s')

# construct cmd line

	CMD="${JAVA_1_8}/java"
	CMD=${CMD}" -jar ${GATK_DIR}/GenomeAnalysisTK.jar"
	CMD=${CMD}" -T CombineVariants"
	CMD=${CMD}" --disable_auto_index_creation_and_locking_when_reading_rods"
	CMD=${CMD}" -R ${GRCH38_REF}"

	# loop to find all the vcf files, index them and add them to the cmd line

		for VCF in \
			$(ls ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.*.FILTERED.GT.REFINED.GRCH38.reject.vcf.gz)
		do
			${TABIX_DIR}/tabix -p vcf ${VCF}

			CMD=${CMD}" --variant ${VCF}"
		done

	CMD=${CMD}" --genotypemergeoption UNSORTED"
	CMD=${CMD}" -o ${CORE_PATH}/${PROJECT_MS}/MULTI_SAMPLE/${PREFIX}.FILTERED.GT.REFINED.LIFTOVER.GRCH38.REJECTED.vcf.gz"

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

END_CAT_VARIANTS=$(date '+%s') # capture time process stops for wall clock tracking purposes.

# write wall clock times to file

	echo ${PROJECT_MS},J01,CAT_REFINED_VARIANTS,${HOSTNAME},${START_CAT_VARIANTS},${END_CAT_VARIANTS} \
	>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}".JOINT.CALL.WALL.CLOCK.TIMES.csv"

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
