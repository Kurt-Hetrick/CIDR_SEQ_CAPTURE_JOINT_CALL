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

START_ALL_INDELS=$(date '+%s') # capture time process starts for wall clock tracking purposes.

# for all samples
# select INDEL,MIXED,MNP,SYMBOLIC sites
# REALLY WE ARE GOING FOR NON-SNP SITES HERE

# construct command line

	CMD="${JAVA_1_8}/java -jar"
	CMD=${CMD}" ${GATK_DIR_4011}/gatk-package-4.0.11.0-local.jar"
	CMD=${CMD}" SelectVariants"
	CMD=${CMD}" --reference ${REF_GENOME}"
	CMD=${CMD}" --variant ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.GT.REFINED.vcf"
	CMD=${CMD}" --output ${CORE_PATH}/${PROJECT_MS}/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/${PREFIX}.GT.REFINED.INDEL.ALL.SAMPLES.vcf"
	CMD=${CMD}" --select-type-to-include INDEL"
	CMD=${CMD}" --select-type-to-include MIXED"
	CMD=${CMD}" --select-type-to-include MNP"
	CMD=${CMD}" --select-type-to-include SYMBOLIC"

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

END_ALL_INDELS=$(date '+%s') # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${PROJECT_MS},J01,ALL_INDELS,${HOSTNAME},${START_ALL_INDELS},${END_ALL_INDELS} \
	>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}.JOINT.CALL.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
