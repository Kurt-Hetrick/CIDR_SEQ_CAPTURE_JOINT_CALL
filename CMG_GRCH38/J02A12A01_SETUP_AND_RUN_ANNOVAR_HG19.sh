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

	PROJECT_SAMPLE=$1
	SM_TAG=$2
	CIDRSEQSUITE_ANNOVAR_JAVA=$3
	CIDRSEQSUITE_DIR_4_0=$4
	CORE_PATH=$5
	CIDRSEQSUITE_PROPS_DIR=$6

# decompress filtered on bait vcf file into a sample temp foler
 
	zcat $CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT_HG19/$SM_TAG"_MS_OnBait.hg19liftover.vcf.gz" \
	>| $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_ANNOVAR_HG19"/$SM_TAG"_MS_OnBait.hg19liftover.vcf"

# run annovar

	$CIDRSEQSUITE_ANNOVAR_JAVA/java -jar \
	-Duser.home=$CIDRSEQSUITE_PROPS_DIR \
	$CIDRSEQSUITE_DIR_4_0/CIDRSeqSuite.jar \
	-pipeline \
	-annovar_directory_annotation \
	$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_ANNOVAR_HG19" \
	$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_ANNOVAR_HG19"

# move the annovar report and dictionary into the sample project /REPORTS/ANNOVAR folder

	du -ah $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_ANNOVAR" \
		| egrep "csv|txt" \
		| awk '{print "mv -v" , $2 , "'$CORE_PATH'" "/" "'$PROJECT_SAMPLE'" "/REPORTS/ANNOVAR_HG19/"}' \
		| bash

# delete the temporary annovar folder

	rm -rvf $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_ANNOVAR_HG19"
