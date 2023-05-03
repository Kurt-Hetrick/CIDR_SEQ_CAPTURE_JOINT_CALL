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
	CIDRSEQSUITE_7_5_0_DIR=$2
	VERACODE_CSV=$3

	CORE_PATH=$4
	PROJECT_SAMPLE=$5
	SM_TAG=$6
	TARGET_BED=$7

START_CONCORDANCE=`date '+%s'`

# decompress the target vcf file into the temporary sub-folder

	zcat $CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_SNV.vcf.gz" \
	>| $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE"/$SM_TAG"_MS_OnTarget_SNV.vcf"

# look for a final report and store it as a variable. if there are multiple ones, then take the newest one

FINAL_REPORT_FILE_TEST=$(ls -tr $CORE_PATH/$PROJECT_SAMPLE/Pretesting/Final_Genotyping_Reports/*$SM_TAG* | tail -n 1)

# if final report exists containing the full sm-tag, then cidrseqsuite magic

if [[ ! -z "$FINAL_REPORT_FILE_TEST" ]]
	then
		FINAL_REPORT=$FINAL_REPORT_FILE_TEST

# if it does not exist, and if the $SM_TAG does not begin with an integer then split $SM_TAG On a @ or -\
# look for a final report that contains that that first element of the $SM_TAG
# bonus feature. if this first tests true but the file still does not exist then cidrseqsuite magic files b/c no file exists

elif [[ $SM_TAG != [0-9]* ]]
	then
	
		HAPMAP=${SM_TAG%[@-]*}

		FINAL_REPORT=$(ls $CORE_PATH/$PROJECT_SAMPLE/Pretesting/Final_Genotyping_Reports/*$HAPMAP* | head -n 1)

else

# both conditions fails then echo the below message and give a dummy value for the $FINAL_REPORT

	echo
	echo At this time, you are looking for a final report that does not exist or fails to meet the current logic for finding a final report.
	echo Please talk to Kurt, because he loves to talk.
	echo

	FINAL_REPORT="FILE_DOES_NOT_EXIST"

fi

# -single_sample_concordance
# Performs concordance between one vcf file and one final report. The vcf must be single sample.
# [1] path_to_vcf_file
# [2] path_to_final_report_file
# [3] path_to_bed_file
# [4] path_to_liftover_file
# [5] path_to_output_directory

CMD=$JAVA_1_8'/java -jar'
CMD=$CMD' '$CIDRSEQSUITE_7_5_0_DIR'/CIDRSeqSuite.jar'
CMD=$CMD' -single_sample_concordance'
CMD=$CMD' '$CORE_PATH'/'$PROJECT_SAMPLE'/TEMP/'$SM_TAG'_MS_CONCORDANCE/'$SM_TAG'_MS_OnTarget_SNV.vcf'
CMD=$CMD' '$FINAL_REPORT
CMD=$CMD' '$TARGET_BED
CMD=$CMD' '$VERACODE_CSV
CMD=$CMD' '$CORE_PATH'/'$PROJECT_SAMPLE'/REPORTS/CONCORDANCE_MS/'

echo $CMD >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
echo $CMD | bash

END_CONCORDANCE=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_SAMPLE",M01,CONCORDANCE,"$HOSTNAME","$START_CONCORDANCE","$END_CONCORDANCE \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".JOINT.CALL.WALL.CLOCK.TIMES.csv"
