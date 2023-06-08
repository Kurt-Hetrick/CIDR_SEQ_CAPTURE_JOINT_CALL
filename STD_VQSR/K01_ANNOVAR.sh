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

	CIDRSEQSUITE_ANNOVAR_JAVA=$1
	CIDRSEQSUITE_DIR_4_0=$2
	CIDRSEQSUITE_PROPS=$3

	CORE_PATH=$4
	PROJECT_MS=$5
	PREFIX=$6

START_ANNOVAR=$(date '+%s') # capture time process starts for wall clock tracking purposes.

# find the vcf and decompress it

	FINAL_MULTI_SAMPLE_VCF=$(ls ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.GT.REFINED{.vcf,.vcf.gz} 2> /dev/null)

	zcat -f ${CORE_PATH}/${PROJECT_MS}/TEMP/${PREFIX}.GT.REFINED{.vcf,.vcf.gz} \
	>| ${CORE_PATH}/${PROJECT_MS}/TEMP/ANNOVAR/${PREFIX}/${PREFIX}.GT.REFINED.vcf

# fyi after -annovar_directory_annotation the arguments are for
# 1. path to vcf file directory
# 2. path to output directory

	CMD="${CIDRSEQSUITE_ANNOVAR_JAVA}/java -jar"
	CMD=${CMD}" -Duser.home=${CIDRSEQSUITE_PROPS}"
	CMD=${CMD}" -Xmx700g"
	CMD=${CMD}" ${CIDRSEQSUITE_DIR_4_0}/CIDRSeqSuite.jar"
	CMD=${CMD}" -pipeline"
	CMD=${CMD}" -annovar_directory_annotation"
	CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/ANNOVAR/${PREFIX}"
	CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/ANNOVAR/${PREFIX}"

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

# move annovar output and dictionary file from temp to REPORTS folder

	du -ah ${CORE_PATH}/${PROJECT_MS}/TEMP/ANNOVAR/${PREFIX} \
		| egrep "csv|txt" \
		| awk '{print "mv -v",$2,"'${CORE_PATH}'" "/" "'${PROJECT_MS}'" "/REPORTS/ANNOVAR/"}' \
		| bash

END_ANNOVAR=$(date '+%s') # capture time process stops for wall clock tracking purposes.

# write out timing metrics to file

	echo ${PROJECT_MS},K01,ANNOVAR,${HOSTNAME},${START_ANNOVAR},${END_ANNOVAR} \
	>> ${CORE_PATH}/${PROJECT_MS}/REPORTS/${PROJECT_MS}.JOINT.CALL.WALL.CLOCK.TIMES.csv

# exit with the signal from the program

	exit ${SCRIPT_STATUS}
