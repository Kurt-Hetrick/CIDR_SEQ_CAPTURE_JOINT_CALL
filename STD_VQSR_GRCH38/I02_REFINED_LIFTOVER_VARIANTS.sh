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

	PICARD_LIFTOVER_CONTAINER=$1

	CORE_PATH=$2
	PROJECT_MS=$3
	PREFIX=$4
	BED_FILE_NAME=$5
	HG19_REF=$6
	HG38_TO_HG19_CHAIN=$7
	SAMPLE_SHEET=$8
		SAMPLE_SHEET_NAME=$(basename ${SAMPLE_SHEET} .csv)
	SUBMIT_STAMP=$9

# liftover from hg38 to hg19

START_LIFTOVER_REFINE_GT=`date '+%s'`

	# construct command line

		CMD="singularity exec ${PICARD_LIFTOVER_CONTAINER} java -jar"
			CMD=${CMD}" /picard/picard.jar"
		CMD=${CMD}" LiftoverVcf"
			CMD=${CMD}" INPUT=${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.${BED_FILE_NAME}.BEDsuperset.VQSR.1KG.ExAC3.REFINED.temp.vcf.gz"
			CMD=${CMD}" REFERENCE_SEQUENCE=${HG19_REF}"
			CMD=${CMD}" CHAIN=${HG38_TO_HG19_CHAIN}"
		CMD=${CMD}" REJECT=${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.${BED_FILE_NAME}.BEDsuperset.VQSR.1KG.ExAC3.REFINED.hg19.reject.vcf.gz"
		CMD=${CMD}" OUTPUT=${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.${BED_FILE_NAME}.BEDsuperset.VQSR.1KG.ExAC3.REFINED.hg19.temp.vcf.gz"
		CMD=${CMD}" &&"
		# remove loci that start with chrUn because cidrseqsuite will crash
		CMD=${CMD}" zegrep -v \"^chrUn\" ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.${BED_FILE_NAME}.BEDsuperset.VQSR.1KG.ExAC3.REFINED.hg19.temp.vcf.gz"
		CMD=${CMD}" >| ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.${BED_FILE_NAME}.BEDsuperset.VQSR.1KG.ExAC3.REFINED.hg19.primary.vcf"

	# write command line to file and execute the command line

		echo ${CMD} >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${PROJECT_MS}_command_lines.txt
		echo >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${PROJECT_MS}_command_lines.txt
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

END_LIFTOVER_REFINE_GT=`date '+%s'`

echo ${PROJECT_MS},J01,REFINE_GT,${HOSTNAME},${START_LIFTOVER_REFINE_GT},${END_LIFTOVER_REFINE_GT} \
>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}.JOINT.CALL.WALL.CLOCK.TIMES.csv

# write out timing metrics to file

	echo ${PROJECT_MS},J01,LIFTOVER_REFINE_GT,${HOSTNAME},${START_LIFTOVER_REFINE_GT},${END_LIFTOVER_REFINE_GT} \
	>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
