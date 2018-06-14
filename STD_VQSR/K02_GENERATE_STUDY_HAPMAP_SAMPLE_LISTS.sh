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

# INPUT VARIABLES

CORE_PATH=$1
PROJECT_MS=$2
PREFIX=$3

# TAKE HEADER IN THE VCF THAT CONTAINS THE SAMPLE NAMES
# CONVERT INTO ROWS
# REMOVE THE STANDARD VCF COLUMNS
# KEEP ONLY THOSE SAMPLE NAMES THAT START WITH "NA" OR "HG"

egrep -m 1 '^#CHROM' $CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf' \
| sed 's/\t/\n/g' \
| awk 'NR>9 {print $0}' \
| egrep '^NA|^HG' \
>| $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/$PREFIX'_hapmap_samples.list'

# TAKE HEADER IN THE VCF THAT CONTAINS THE SAMPLE NAMES
# CONVERT INTO ROWS
# REMOVE THE STANDARD VCF COLUMNS
# KEEP ONLY THOSE SAMPLE NAMES THAT DO NOT START WITH "NA" OR "HG"

egrep -m 1 '^#CHROM' $CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf' \
| sed 's/\t/\n/g' \
| awk 'NR>9 {print $0}' \
| egrep -v '^NA|^HG' \
>| $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/$PREFIX'_study_samples.list'
