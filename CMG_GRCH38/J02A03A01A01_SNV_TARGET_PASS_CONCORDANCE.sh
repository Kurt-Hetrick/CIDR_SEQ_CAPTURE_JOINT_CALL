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
	CIDRSEQSUITE_7_5_0_DIR=$2
	VERACODE_CSV=$3
	
	CORE_PATH=$4
	PROJECT_SAMPLE=$5
	SM_TAG=$6
	TARGET_BED=$7
		TARGET_BED_NAME=(`basename $TARGET_BED .bed`)
	REF_GENOME=$8
		REF_DIR=$(dirname $REF_GENOME)
		REF_BASENAME=$(basename $REF_GENOME | sed 's/.fasta//g ; s/.fa//g')
	PICARD_DIR=$9
	H19_DICT=${10}
	HG38_TO_HG19_CHAIN=${11}
	BEDTOOLS_DIR=${12}

# MAKE PICARD INTERVAL FILES (1-based start)

	(grep "^@SQ" $REF_DIR/$REF_BASENAME".dict" \
		; awk 1 $TARGET_BED \
			| sed -r 's/\r//g ; s/[[:space:]]+/\t/g' \
			| awk 'BEGIN {OFS="\t"} {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}') \
	>| $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG".OnTarget.picard.bed"

# LIFTOVER PICARD TARGET INTERVAL LIST BACK TO HG19
# this is b/c we suck and are stuck on RHEL6 and current version of liftover needs GLIBC 2.14+

	$JAVA_1_8/java -jar \
	$PICARD_DIR/picard.jar \
	LiftOverIntervalList \
	INPUT=$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG".OnTarget.picard.bed" \
	OUTPUT=$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG".OnTarget.hg19.interval_list" \
	SEQUENCE_DICTIONARY=$H19_DICT \
	CHAIN=$HG38_TO_HG19_CHAIN

		# CONVERT HG19 CONVERT INTERVAL LIST BACK TO A BED FILE
				
			grep -v "^@" $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG".OnTarget.hg19.interval_list" \
			| awk 'BEGIN {OFS="\t"} {print $1,($2-1),$3}' \
			| $BEDTOOLS_DIR/bedtools merge -i - \
			>| $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"-"TARGET_BED_NAME".lift.hg19.bed"

# decompress the target vcf file into the temporary sub-folder

	zegrep -v "^chrUn" $CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET_HG19/$SM_TAG"_MS_OnTarget_SNV.hg19liftover.vcf.gz" \
	>| $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE"/$SM_TAG"_MS_OnTarget_SNV.hg19liftover.vcf"

# look for a final report and store it as a variable. if there are multiple ones, then take the newest one

	FINAL_REPORT_FILE_TEST=$(ls -tr $CORE_PATH/$PROJECT_SAMPLE/Pretesting/Final_Genotyping_Reports/*$SM_TAG* | tail -n 1)

# if final report exists containing the full sm-tag, then cidrseqsuite magic

if [[ ! -z "$FINAL_REPORT_FILE_TEST" ]];then
	
		FINAL_REPORT=$FINAL_REPORT_FILE_TEST

# if it does not exist, and if the $SM_TAG does not begin with an integer then split $SM_TAG On a @ or -\
# look for a final report that contains that that first element of the $SM_TAG
# bonus feature. if this first tests true but the file still does not exist then cidrseqsuite magic files b/c no file exists

	elif [[ $SM_TAG != [0-9]* ]]; then
		
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

# I should fix the target bed file on the fly here...but I'm not at the moment.

	$JAVA_1_8/java -jar \
	$CIDRSEQSUITE_7_5_0_DIR/CIDRSeqSuite.jar \
	-single_sample_concordance \
	$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE"/$SM_TAG"_MS_OnTarget_SNV.hg19liftover.vcf" \
	$FINAL_REPORT \
	$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"-"TARGET_BED_NAME".lift.hg19.bed" \
	$VERACODE_CSV \
	$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE"

echo \
$JAVA_1_8/java -jar \
$CIDRSEQSUITE_7_5_0_DIR/CIDRSeqSuite.jar \
-single_sample_concordance \
$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE"/$SM_TAG"_MS_OnTarget_SNV.hg19liftover.vcf" \
$FINAL_REPORT \
$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"-"TARGET_BED_NAME".lift.hg19.bed" \
$VERACODE_CSV \
$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE" \
>> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

	mv -v $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE"/$SM_TAG"_concordance.csv" \
	$CORE_PATH/$PROJECT_SAMPLE/REPORTS/CONCORDANCE_MS/$SM_TAG"_concordance.csv"

	mv -v $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE"/$SM_TAG"_discordant_calls.txt" \
	$CORE_PATH/$PROJECT_SAMPLE/REPORTS/CONCORDANCE_MS/$SM_TAG"_discordant_calls.txt"
