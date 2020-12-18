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
# redirecting stderr/stdout to file as a log.

	set

	echo

# INPUT PARAMETERS

	JAVA_1_8=$1
	LAB_QC_DIR=$2
	CORE_PATH=$3

	PROJECT_MS=$4
	SAMPLE_SHEET=$5
		SAMPLE_SHEET_NAME=`basename $SAMPLE_SHEET .csv`

START_LAB_PREP_METRICS=`date '+%s'`

# Generates a QC report for lab specific metrics including Physique Report, Samples Table, Sequencer XML data, Pca and Phoenix. Does not check if samples are dropped.
	# [1] path_to_sample_sheet
	# [2] path_to_seq_proj ($CORE_PATH)
	# [3] path_to_output_file

		$JAVA_1_8/java -jar $LAB_QC_DIR/EnhancedSequencingQCReport.jar \
		-lab_qc_metrics \
		$SAMPLE_SHEET \
		$CORE_PATH \
		$CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".LAB_PREP_METRICS.csv"

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

	# if exit does not equal 0 then exit with whatever the exit signal is at the end.
	# also write to file that this job failed

		if [ "$SCRIPT_STATUS" -ne 0 ]
		 then
			echo $PROJECT $HOSTNAME $JOB_NAME $USER $SCRIPT_STATUS $SGE_STDERR_PATH \
			>> $CORE_PATH/$PROJECT/TEMP/$SAMPLE_SHEET_NAME"_"$SUBMIT_STAMP"_ERRORS.txt"
			exit $SCRIPT_STATUS
		fi

END_LAB_PREP_METRICS=`date '+s'`

(head -n 1 $CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".LAB_PREP_METRICS.csv" \
	| awk '{print $0 ",EPOCH_TIME"}' ; \
	awk 'NR>1 {print $0 "," "'$START_LAB_PREP_METRICS'"}' \
	$CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".LAB_PREP_METRICS.csv" \
	| sort -k 1,1 ) \
	>| $CORE_PATH/$PROJECT_MS/REPORTS/LAB_PREP_REPORTS_MS/$SAMPLE_SHEET_NAME".LAB_PREP_METRICS.csv"

echo $PROJECT_MS,X.01,LAB_QC_PREP_METRICS,$HOSTNAME,$START_LAB_PREP_METRICS_METRICS,$END_LAB_PREP_METRICS \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $LAB_QC_DIR/EnhancedSequencingQCReport.jar \
-lab_qc_metrics \
$CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME"_"$START_LAB_PREP_METRICS".csv" \
$CORE_PATH \
$CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".LAB_PREP_METRICS.csv" \
>> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"

echo $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"

# exit with the signal from the program

	exit $SCRIPT_STATUS
