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

	TABIX_DIR=$1

	CORE_PATH=$2
	PROJECT_SAMPLE=$3
	PROJECT_MS=$4
	SM_TAG=$5
	PREFIX=$6

START_SELECT_ALL_SAMPLE=`date '+%s'`

SAMPLE_COLUMN=(`egrep -m 1 ^#CHROM $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.vcf" \
	| sed 's/\t/\n/g' \
	| cat -n \
	| grep -w $SM_TAG \
	| sed 's/^ *//g' \
	| cut -f 1`)

# Extract out sample, remove non-passing, non-variant

	awk \
	-v SAMPLE_COLUMN="$SAMPLE_COLUMN" \
	'BEGIN {FS="\t";OFS="\t"} \
	$SAMPLE_COLUMN!~/^[0]\/[0]/&&$7!="LowQual"&&$SAMPLE_COLUMN!~/^[.]\/[.]/&&$7!~"VQSR" \
	{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$SAMPLE_COLUMN}' \
	$CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.SNP.HardFiltered.INDEL.1KG.ExAC3.REFINED.vcf" \
	>| $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_OnBait.vcf"

# compress with bgzip

	$TABIX_DIR/bgzip -c $CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_OnBait.vcf" \
	>| $CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz"

# index with tabix

	$TABIX_DIR/tabix -p vcf \
	$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz"

END_SELECT_ALL_SAMPLE=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT_SAMPLE",K01,SELECT_ALL_SAMPLE,"$HOSTNAME","$START_SELECT_ALL_SAMPLE","$END_SELECT_ALL_SAMPLE \
>> $CORE_PATH/$PROJECT_SAMPLE/REPORTS/$PROJECT_SAMPLE".JOINT.CALL.WALL.CLOCK.TIMES.csv"

ls $CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz.tbi"
