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

# create a blank lane b/w the output variables and the program logging output

	echo

# INPUT VARIABLES

	JAVA_1_8=$1
	PICARD_DIR=$2
	
	CORE_PATH=$3
	PROJECT_MS=$4
	PREFIX=$5
	GRCH38_REF=$6
	HG19_TO_GRCH38_CHAIN=$7

# liftover from hg19 to GRCh38

START_LIFTOVER_MS_HG38=$(date '+%s') # capture time process starts for wall clock tracking purposes.

# construct cmd line

	CMD="${JAVA_1_8}/java -jar"
	CMD=${CMD}" ${PICARD_DIR}/picard.jar"
	CMD=${CMD}" LiftoverVcf"
	CMD=${CMD}" INPUT=${CORE_PATH}/${PROJECT_MS}/MULTI_SAMPLE/${PREFIX}.FILTERED.HG19.LIFTOVER.vcf.gz"
	CMD=${CMD}" OUTPUT=${CORE_PATH}/${PROJECT_MS}/MULTI_SAMPLE/${PREFIX}.FILTERED.GRCh38.LIFTOVER.vcf.gz"
	CMD=${CMD}" REJECT=${CORE_PATH}/${PROJECT_MS}/MULTI_SAMPLE/${PREFIX}.FILTERED.GRCh38.LIFTOVER.REJECTED.vcf.gz"
	CMD=${CMD}" REFERENCE_SEQUENCE=${GRCH38_REF}"
	CMD=${CMD}" CHAIN=${HG19_TO_GRCH38_CHAIN}"

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

END_LIFTOVER_MS_HG38=$(date '+%s') # capture time process stops for wall clock tracking purposes.

# write wall clock times to file

	echo $SM_TAG_${PROJECT_MS},I.01,LIFTOVER_INITIAL_VCF_HG38,${HOSTNAME},${START_LIFTOVER_MS_HG38},${END_LIFTOVER_MS_HG38} \
	>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
