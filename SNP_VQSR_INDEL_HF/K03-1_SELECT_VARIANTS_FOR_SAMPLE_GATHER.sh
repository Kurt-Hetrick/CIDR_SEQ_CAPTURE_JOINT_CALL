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
	GATK_DIR=$2
	SAMPLE_REF_GENOME=$3

	CORE_PATH=$4
	PROJECT_MS=$5
	PROJECT_SAMPLE=$6
	SM_TAG=$7
	BAIT_BED=$8

## -----CONCATENATE SCATTERED g.vcf FILES INTO A SINGLE GRCh37 reference sorted g.vcf file-----

# Start with creating a *list file, reference sorted, to put into --variant.
# Assumption is that this is a correctly sorted GRCh37 reference file as the input reference used

# Put the autosome into a file, sort numerically

sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
| sed -r 's/[[:space:]]+/\t/g' \
| cut -f 1 \
| sort \
| uniq \
| awk '$1~/^[0-9]/' \
| sort -k1,1n \
| awk '{print "'$CORE_PATH'" "/" "'$PROJECT_MS'" "/TEMP/" "'$SM_TAG'" "_SCATTER/"  "'$SM_TAG'" "_MS_OnBait." $1 ".vcf.gz"}' \
>| $CORE_PATH/$PROJECT_MS/TEMP/$SM_TAG".FINAL.list"

# Append X if present

sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
| sed -r 's/[[:space:]]+/\t/g' \
| cut -f 1 \
| sort \
| uniq \
| awk '$1=="X"' \
| awk '{print "'$CORE_PATH'" "/" "'$PROJECT_MS'" "/TEMP/" "'$SM_TAG'" "_SCATTER/"  "'$SM_TAG'" "_MS_OnBait." $1 ".vcf.gz"}' \
>> $CORE_PATH/$PROJECT_MS/TEMP/$SM_TAG".FINAL.list"

# Append Y if present

sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
| sed -r 's/[[:space:]]+/\t/g' \
| cut -f 1 \
| sort \
| uniq \
| awk '$1=="Y"' \
| awk '{print "'$CORE_PATH'" "/" "'$PROJECT_MS'" "/TEMP/" "'$SM_TAG'" "_SCATTER/"  "'$SM_TAG'" "_MS_OnBait." $1 ".vcf.gz"}' \
>> $CORE_PATH/$PROJECT_MS/TEMP/$SM_TAG".FINAL.list"

# Append MT if present

sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
| sed -r 's/[[:space:]]+/\t/g' \
| cut -f 1 \
| sort \
| uniq \
| awk '$1=="MT"' \
| awk '{print "'$CORE_PATH'" "/" "'$PROJECT_MS'" "/TEMP/" "'$SM_TAG'" "_SCATTER/"  "'$SM_TAG'" "_MS_OnBait." $1 ".vcf.gz"}' \
>> $CORE_PATH/$PROJECT_MS/TEMP/$SM_TAG".FINAL.list"

## Gather the per chromosome FINAL files

START_FINAL_SAMPLE_VCF_GATHER=`date '+%s'`

$JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $SAMPLE_REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT_MS/TEMP/$SM_TAG".FINAL.list" \
--outputFile $CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz"

END_FINAL_SAMPLE_VCF_GATHER=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",K03-1,FINAL_SAMPLE_VCF_GATHER,"$HOSTNAME","$START_FINAL_SAMPLE_VCF_GATHER","$END_FINAL_SAMPLE_VCF_GATHER \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS"JOINT.CALL.WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -cp $GATK_DIR/GenomeAnalysisTK.jar \
org.broadinstitute.gatk.tools.CatVariants \
-R $REF_GENOME \
--assumeSorted \
--variant $CORE_PATH/$PROJECT_MS/TEMP/$SM_TAG".FINAL.list" \
--outputFile $CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz" \
>> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT_SAMPLE/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

# if file is not present exit !=0

ls $CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait.vcf.gz.tbi"
