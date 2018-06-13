#! /bin/bash

PROJECT_MS=$1 # the project where the multi-sample vcf is being written to
SAMPLE_SHEET=$2 # full/relative path to the sample sheet
PREFIX=$3 # prefix name that you want to give the multi-sample vcf
NUMBER_OF_BED_FILES=$4 # scatter count, if not supplied then the default is what is below.

# if there is no 4 the argument present then use the number for the scatter count
if [[ ! $NUMBER_OF_BED_FILES ]]
	then
	NUMBER_OF_BED_FILES=500
fi

# datamash is used in the submitter
# gcc is so that it can be pushed out to the compute nodes via qsub (-V)
module load datamash
module load gcc/5.1.0

# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED
SCRIPT_DIR="/mnt/research/tools/LINUX/00_GIT_REPO_KURT/CIDR_SEQ_CAPTURE_JOINT_CALL/STD_VQSR"

# Directory where sequencing projects are located
CORE_PATH="/mnt/research/active"

# Generate a list of active queue and remove the ones that I don't want to use
QUEUE_LIST=`qstat -f -s r \
 | egrep -v "^[0-9]|^-|^queue" \
 | cut -d @ -f 1 \
 | sort \
 | uniq \
 | egrep -v "bigmem.q|all.q|cgc.q|programmers.q|rhel7.q" \
 | datamash collapse 1 \
 | awk '{print $1}'`

# EVENTUALLY I WANT THIS SET UP AS AN OPTION WITH A DEFAULT OF X

PRIORITY="-15"

# eventually, i want to push this out to something...maybe in the vcf file header.
PIPELINE_VERSION=`git --git-dir=$SCRIPT_DIR/../.git --work-tree=$SCRIPT_DIR/.. log --pretty=format:'%h' -n 1`

#####################
# PIPELINE PROGRAMS #
#####################

JAVA_1_8="/mnt/linuxtools/JAVA/jdk1.8.0_73/bin"

BEDTOOLS_DIR="/mnt/linuxtools/BEDTOOLS/bedtools-2.22.0/bin"
GATK_DIR="/mnt/linuxtools/GATK/GenomeAnalysisTK-3.7"
SAMTOOLS_0118_DIR="/mnt/linuxtools/SAMTOOLS/samtools-0.1.18"
	# Becasue I didn't want to go through compiling this yet for version 1.6...I'm hoping that Keith will eventually do a full OS install of RHEL7 instead of his
	# typical stripped down installations so I don't have to install multiple libraries again
TABIX_DIR="/mnt/linuxtools/TABIX/tabix-0.2.6"
CIDR_SEQSUITE_JAVA_DIR="/mnt/linuxtools/JAVA/jre1.7.0_45/bin"
CIDR_SEQSUITE_6_1_1_DIR="/mnt/linuxtools/CIDRSEQSUITE/6.1.1"

##################
# PIPELINE FILES #
##################

HAPMAP_VCF="/mnt/research/tools/PIPELINE_FILES/GATK_resource_bundle/2.8/b37/hapmap_3.3.b37.vcf"
OMNI_VCF="/mnt/research/tools/PIPELINE_FILES/GATK_resource_bundle/2.8/b37/1000G_omni2.5.b37.vcf"
ONEKG_SNPS_VCF="/mnt/research/tools/PIPELINE_FILES/GATK_resource_bundle/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf"
DBSNP_138_VCF="/mnt/research/tools/PIPELINE_FILES/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.vcf"
ONEKG_INDELS_VCF="/mnt/research/tools/PIPELINE_FILES/GATK_resource_bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
P3_1KG="/mnt/research/tools/PIPELINE_FILES/GRCh37_aux_files/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
ExAC="/mnt/research/tools/PIPELINE_FILES/GRCh37_aux_files/ExAC.r0.3.sites.vep.vcf.gz"
KNOWN_SNPS="/mnt/research/tools/PIPELINE_FILES/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf"
VERACODE_CSV="/mnt/linuxtools/CIDRSEQSUITE/Veracode_hg18_hg19.csv"

##########################
##### PROJECT SET-UP #####
##########################

## This checks to see if bed file directory has been created from a previous run.  If so, remove it to not interfere with current run
if [ -d $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT ]
then
	rm -rf $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT
fi

# MAKE THE FOLLOWING FOLDERS IN THE PROJECT WHERE THE MULTI-SAMPLE VCF IS GOING TO
mkdir -p $CORE_PATH/$PROJECT_MS/{LOGS,COMMAND_LINES}
mkdir -p $CORE_PATH/$PROJECT_MS/TEMP/{BED_FILE_SPLIT,AGGREGATE}
mkdir -p $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/
mkdir -p $CORE_PATH/$PROJECT_MS/GVCF/AGGREGATE

# grab the reference genome file, dbsnp file and bait bed file for the "project"
## !should do a check here to make sure that there is only one record...!

CREATE_PROJECT_INFO_ARRAY ()
{
PROJECT_INFO_ARRAY=(`sed 's/\r//g' $SAMPLE_SHEET \
	| awk 'BEGIN{FS=","} NR>1 {print $12,$18,$16}' \
	| sed 's/,/\t/g' \
	| sort -k 1,1 \
	| awk '{print "'$PROJECT_MS'",$1,$2,$3}' \
	| sort \
	| uniq`)

SEQ_PROJECT=${PROJECT_INFO_ARRAY[0]} # same as $PROJECT_MS...
REF_GENOME=${PROJECT_INFO_ARRAY[1]} # field 12 from the sample sheet
PROJECT_DBSNP=${PROJECT_INFO_ARRAY[2]} # field 18 from the sample sheet
PROJECT_BAIT_BED=${PROJECT_INFO_ARRAY[3]} # field 16 from the sample sheet
}

CREATE_GVCF_LIST()
{
# count how many unique sample id's (with project) are in the sample sheet.
TOTAL_SAMPLES=(`awk 'BEGIN{FS=","} NR>1{print $1,$8}' $SAMPLE_SHEET \
	| sort \
	| uniq \
	| wc -l`)

# find all of the gvcf files write all of the full paths to a *samples.gvcf.list file.
awk 'BEGIN{FS=","} NR>1{print $1,$8}' $SAMPLE_SHEET \
 | sort \
 | uniq \
 | awk 'BEGIN{OFS="/"}{print "ls " "'$CORE_PATH'",$1,"GVCF",$2".g.vcf*"}' \
 | bash \
 | egrep -v "idx|tbi" \
>| $CORE_PATH'/'$PROJECT_MS'/'$TOTAL_SAMPLES'.samples.gvcf.list'

# STORE THE GVCF LIST FILE PATH AS A VARIABLE
GVCF_LIST=(`echo $CORE_PATH'/'$PROJECT_MS'/'$TOTAL_SAMPLES'.samples.gvcf.list'`)
}

#################

FORMAT_AND_SCATTER_BAIT_BED() 
{
BED_FILE_PREFIX=(`echo SPLITTED_BED_FILE_`)

awk 1 $PROJECT_BAIT_BED | sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
>| $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed

(awk '$1~/^[0-9]/' $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed | sort -k1,1n -k2,2n ; \
 awk '$1=="X"' $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed | sort -k 2,2n ; \
 awk '$1=="Y"' $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed | sort -k 2,2n ; \
 awk '$1=="MT"' $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed | sort -k 2,2n) \
>| $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_AND_SORTED_BED_FILE.bed

# Determining how many records will be in each mini-bed file.
# The +1 at the end is to round up the number of records per mini-bed file to ensure all records are captured.
# So the last mini-bed file will be smaller.
## IIRC. this statement isn't really true, but I don't feel like figuring it out right now. KNH
INTERVALS_DIVIDED=`wc -l $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_AND_SORTED_BED_FILE.bed \
	| awk '{print $1"/""'$NUMBER_OF_BED_FILES'"}' \
	| bc \
	| awk '{print $0+1}'`

split -l $INTERVALS_DIVIDED -a 4 -d  $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_AND_SORTED_BED_FILE.bed \
$CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/$BED_FILE_PREFIX

ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/$BED_FILE_PREFIX* | awk '{print "mv",$0,$0".bed"}' | bash
}

CREATE_PROJECT_INFO_ARRAY
MAKE_PROJ_DIR_TREE
FORMAT_AND_SCATTER_BAIT_BED
CREATE_GVCF_LIST

###############################################################################################################
###############################################################################################################
###############################################################################################################

############################################################################
#################Start of Combine Gvcf Functions############################
############################################################################

COMBINE_GVCF()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N 'A01_COMBINE_GVCF_'$PROJECT_MS'_'$BED_FILE_NAME \
		-o $CORE_PATH/$PROJECT_MS/LOGS/A01_COMBINE_GVCF_$BED_FILE_NAME.log \
	$SCRIPT_DIR/A01_COMBINE_GVCF.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$GVCF_LIST \
		$PREFIX \
		$BED_FILE_NAME
}

GENOTYPE_GVCF()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N B02_GENOTYPE_GVCF_$PROJECT_MS'_'$BED_FILE_NAME \
		-o $CORE_PATH/$PROJECT_MS/LOGS/B02_GENOTYPE_GVCF_$BED_FILE_NAME.log \
		-hold_jid A01_COMBINE_GVCF_$PROJECT_MS'_'$BED_FILE_NAME \
	$SCRIPT_DIR/B02_GENOTYPE_GVCF.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX \
		$BED_FILE_NAME
}

VARIANT_ANNOTATOR()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N C03_VARIANT_ANNOTATOR_$PROJECT_MS'_'$BED_FILE_NAME \
		-o $CORE_PATH/$PROJECT_MS/LOGS/C03_VARIANT_ANNOTATOR_$BED_FILE_NAME.log \
		-hold_jid B02_GENOTYPE_GVCF_$PROJECT_MS'_'$BED_FILE_NAME \
	$SCRIPT_DIR/C03_VARIANT_ANNOTATOR.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX \
		$BED_FILE_NAME \
		$PROJECT_DBSNP
}

GENERATE_CAT_VARIANTS_HOLD_ID()
{
	CAT_VARIANTS_HOLD_ID=$CAT_VARIANTS_HOLD_ID'C03_VARIANT_ANNOTATOR_'$PROJECT_MS'_'$BED_FILE_NAME','
}

for BED_FILE in $(ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/SPLITTED_BED_FILE*);
do
	BED_FILE_NAME=$(basename $BED_FILE .bed)
	COMBINE_GVCF
	GENOTYPE_GVCF
	VARIANT_ANNOTATOR
	GENERATE_CAT_VARIANTS_HOLD_ID
done

##############################################################################
#####################End of Combine Gvcf Functions############################
##############################################################################

##############################################################################
##################Start of VQSR and Refinement Functions######################
##############################################################################

CAT_VARIANTS()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N D04_CAT_VARIANTS_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/D04_CAT_VARIANTS.log \
		-hold_jid $CAT_VARIANTS_HOLD_ID \
	$SCRIPT_DIR/D04_CAT_VARIANTS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

VARIANT_RECALIBRATOR_SNV()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \ 
	-N E05A_VARIANT_RECALIBRATOR_SNV_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/E05A_VARIANT_RECALIBRATOR_SNV.log \
		-hold_jid D04_CAT_VARIANTS_$PROJECT_MS \
	$SCRIPT_DIR/E05A_VARIANT_RECALIBRATOR_SNV.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$HAPMAP_VCF \
		$OMNI_VCF \
		$ONEKG_SNPS_VCF \
		$DBSNP_138_VCF \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

VARIANT_RECALIBRATOR_INDEL()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \ 
	-N E05B_VARIANT_RECALIBRATOR_INDEL_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/E05B_VARIANT_RECALIBRATOR_INDEL.log \
		-hold_jid D04_CAT_VARIANTS_$PROJECT_MS \
	$SCRIPT_DIR/E05B_VARIANT_RECALIBRATOR_INDEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$ONEKG_INDELS_VCF \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

APPLY_RECALIBRATION_SNV()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N F06_APPLY_RECALIBRATION_SNV_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/F06_APPLY_RECALIBRATION_SNV.log \
		-hold_jid E05A_VARIANT_RECALIBRATOR_SNV_$PROJECT_MS \
	$SCRIPT_DIR/F06_APPLY_RECALIBRATION_SNV.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

APPLY_RECALIBRATION_INDEL()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \ 
	-N G07_APPLY_RECALIBRATION_INDEL_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/G07_APPLY_RECALIBRATION_INDEL.log \
		-hold_jid F06_APPLY_RECALIBRATION_SNV_$PROJECT_MS,\
 E05B_VARIANT_RECALIBRATOR_INDEL_$PROJECT_MS \
	$SCRIPT_DIR/G07_APPLY_RECALIBRATION_INDEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

# NO LONGER NEEDED.

# BGZIP_AND_TABIX_RECAL_VCF()
# {
# 	echo \
# 	qsub \
# 		-S /bin/bash \
# 		-cwd \
# 		-V \
# 		-q $QUEUE_LIST \
# 		-p $PRIORITY \
# 		-j y \ 
# 	-N H08A_BGZIP_AND_TABIX_RECAL_VCF_$PROJECT_MS \
# 		-o $CORE_PATH/$PROJECT_MS/LOGS/H08A_BGZIP_AND_TABIX_RECAL_VCF.log \
# 		-hold_jid G07_APPLY_RECALIBRATION_INDEL_$PROJECT_MS,\
#  F06_APPLY_RECALIBRATION_SNV_$PROJECT_MS \
# 	$SCRIPT_DIR/H08A_BGZIP_AND_TABIX_RECAL_VCF.sh \
# 		$TABIX_DIR \
# 		$CORE_PATH \
# 		$PROJECT_MS \
# 		$PREFIX
# }


# Change this into a scatter gather.
CALCULATE_GENOTYPE_POSTERIORS()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N H08B_CALCULATE_GENOTYPE_POSTERIORS_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/H08B_CALCULATE_GENOTYPE_POSTERIORS.log \
		-hold_jid G07_APPLY_RECALIBRATION_INDEL_$PROJECT_MS \
	$SCRIPT_DIR/H08B_CALCULATE_GENOTYPE_POSTERIORS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$P3_1KG \
		$ExAC \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

# this would remain in the scatter
VARIANT_ANNOTATOR_REFINED()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N I09_VARIANT_ANNOTATOR_REFINED_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/I09_VARIANT_ANNOTATOR_REFINED.log \
		-hold_jid H08B_CALCULATE_GENOTYPE_POSTERIORS_$PROJECT_MS \
	$SCRIPT_DIR/I09_VARIANT_ANNOTATOR_REFINED.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$PROJECT_DBSNP \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

# there would be a gather in here

# THIS STEP WILL BE REDUNDANT

BGZIP_AND_TABIX_REFINED_VCF()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N J10_BGZIP_AND_TABIX_REFINED_VCF_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/J10_BGZIP_AND_TABIX_REFINED_VCF.log \
		-hold_jid I09_VARIANT_ANNOTATOR_REFINED_$PROJECT_MS \
	$SCRIPT_DIR/J10_BGZIP_AND_TABIX_REFINED_VCF.sh \
		$TABIX_DIR \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

CAT_VARIANTS
VARIANT_RECALIBRATOR_SNV
VARIANT_RECALIBRATOR_INDEL
APPLY_RECALIBRATION_SNV
APPLY_RECALIBRATION_INDEL
# BGZIP_AND_TABIX_RECAL_VCF
CALCULATE_GENOTYPE_POSTERIORS # SHOULD CONVERT TO A SCATTER GATHER
VARIANT_ANNOTATOR_REFINED # THIS SHOULD REMAIN SCATTERED
# THERE SHOULD BE A GATHER STEP
# PLACE TO ADD ANNOVAR
BGZIP_AND_TABIX_REFINED_VCF # THIS SHOULD BE REDUNDANT

###########################################################################
#################End of VQSR and Refinement Functions######################
###########################################################################
#
###########################################################################
###################Start of Vcf Splitter Functions#########################
###########################################################################

CREATE_SAMPLE_INFO_ARRAY ()
{
	SAMPLE_INFO_ARRAY=(`sed 's/\r//g' $SAMPLE_SHEET \
		| awk 'BEGIN{FS=","} NR>1 {print $1,$8,$17,$15,$18,$12}' \
		| sed 's/,/\t/g' \
		| sort -k 8,8 \
		| uniq \
		| awk '$2=="'$SAMPLE'" {print $1,$2,$3,$4,$5,$6}'`)
	
	PROJECT_SAMPLE=${SAMPLE_INFO_ARRAY[0]}
	SM_TAG=${SAMPLE_INFO_ARRAY[1]}
	TARGET_BED=${SAMPLE_INFO_ARRAY[2]}
	TITV_BED=${SAMPLE_INFO_ARRAY[3]}
	DBSNP=${SAMPLE_INFO_ARRAY[4]} #Not used unless we implement HC_BAM
	SAMPLE_REF_GENOME=${SAMPLE_INFO_ARRAY[5]}
	
	UNIQUE_ID_SM_TAG=$(echo $SM_TAG | sed 's/@/_/g') # If there is an @ in the qsub or holdId name it breaks
}

# DOES THIS NEEDS GO HERE?
# do need a temp/sm_tag folder

MAKE_PROJ_DIR_TREE ()
{
	mkdir -p \
	$CORE_PATH/$SEQ_PROJECT/{TEMP,LOGS,COMMAND_LINES} \
	$CORE_PATH/$SEQ_PROJECT/INDEL/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
	$CORE_PATH/$SEQ_PROJECT/SNV/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
	$CORE_PATH/$SEQ_PROJECT/MIXED/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
	$CORE_PATH/$SEQ_PROJECT/VCF/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
	$CORE_PATH/$SEQ_PROJECT/REPORTS/{ANNOVAR,QC_REPORT_PREP_MS,QC_REPORTS,LAB_PREP_REPORTS,TI_TV_MS,CONCORDANCE_MS} \
	$CORE_PATH/$SEQ_PROJECT/TEMP/$SM_TAG
}


# NEED TO LOOK AT HOLD ID
SELECT_PASSING_VARIANTS_PER_SAMPLE()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_MS/LOGS/K11_SELECT_VARIANTS_FOR_SAMPLE_$SM_TAG.log \
		-hold_jid I09_VARIANT_ANNOTATOR_REFINED_$PROJECT \
	$SCRIPT_DIR/K11_SELECT_VARIANTS_FOR_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$PREFIX
}

# THIS IS UNNECESSARY
# REMOVE

# BGZIP_AND_TABIX_SAMPLE_VCF()
# {
# 	echo \
# 	qsub \
# 		-S /bin/bash \
# 		-cwd \
# 		-V \
# 		-q $QUEUE_LIST \
# 		-p $PRIORITY \
# 		-j y \ 
# 	-N L12A_BGZIP_AND_TABIX_SAMPLE_VCF_$UNIQUE_ID_SM_TAG \
# 		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
# 		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12A_BGZIP_AND_TABIX_SAMPLE_VCF_$SAMPLE.log \
# 	$SCRIPT_DIR/L12A_BGZIP_AND_TABIX_SAMPLE_VCF.sh \
# 		$TABIX_DIR \
# 		$CORE_PATH \
# 		$PROJECT_SAMPLE \
# 		$SM_TAG
# }

PASSING_VARIANTS_ON_TARGET_BY_SAMPLE()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12B_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12B_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
	$SCRIPT_DIR/L12B_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TARGET_BED
}

PASSING_SNVS_ON_BAIT_BY_SAMPLE()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12C_PASSING_SNVS_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12C_PASSING_SNVS_ON_BAIT_BY_SAMPLE_$SAMPLE.log \
		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12C_PASSING_SNVS_ON_BAIT_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG
}

PASSING_SNVS_ON_TARGET_BY_SAMPLE()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12D_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12D_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12D_PASSING_SNVS_ON_TARGET_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TARGET_BED
}

PASSING_INDELS_ON_BAIT_BY_SAMPLE()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12E_PASSING_INDELS_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12E_PASSING_INDELS_ON_BAIT_BY_SAMPLE_$SAMPLE.log \
		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12E_PASSING_INDELS_ON_BAIT_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG
}

PASSING_INDELS_ON_TARGET_BY_SAMPLE()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12F_PASSING_INDELS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12F_PASSING_INDELS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12F_PASSING_INDELS_ON_TARGET_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TARGET_BED
}

PASSING_SNVS_TITV_ALL()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12G_PASSING_SNVS_TITV_ALL_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12G_PASSING_SNVS_TITV_ALL_$SAMPLE.log \
		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12G_PASSING_SNVS_TITV_ALL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TITV_BED
}

TITV_ALL()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12G-1_TITV_ALL_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12G-1_TITV_ALL_$SAMPLE.log \
		-hold_jid L12G_PASSING_SNVS_TITV_ALL_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12G-1_TITV_ALL.sh \
		$SAMTOOLS_0118_DIR \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG
}

PASSING_SNVS_TITV_KNOWN()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12H_PASSING_SNVS_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12H_PASSING_SNVS_TITV_KNOWN_$SAMPLE.log \
		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12H_PASSING_SNVS_TITV_KNOWN.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$KNOWN_SNPS \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TITV_BED
}

TITV_KNOWN()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12H-1_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12H-1_TITV_KNOWN_$SAMPLE.log \
		-hold_jid L12H_PASSING_SNVS_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12H-1_TITV_KNOWN.sh \
		$SAMTOOLS_0118_DIR \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG
}

PASSING_SNVS_TITV_NOVEL()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12I_PASSING_SNVS_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12I_PASSING_SNVS_TITV_NOVEL_$SAMPLE.log \
		-hold_jid K11_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12I_PASSING_SNVS_TITV_NOVEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$KNOWN_SNPS \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TITV_BED
}

TITV_NOVEL()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12I-1_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12I-1_TITV_NOVEL_$SAMPLE.log \
		-hold_jid L12I_PASSING_SNVS_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12I-1_TITV_NOVEL.sh \
		$SAMTOOLS_0118_DIR \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG
}

# reconfigure this to match up the current qc pipeline

CONCORDANCE_ON_TARGET_PER_SAMPLE()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N L12D-1_CONCORDANCE_ON_TARGET_PER_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/L12D-1_CONCORDANCE_ON_TARGET_PER_SAMPLE_$SAMPLE.log \
		-hold_jid L12D_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/L12D-1_CONCORDANCE_ON_TARGET_PER_SAMPLE.sh \
		$CIDR_SEQSUITE_JAVA_DIR \
		$CIDR_SEQSUITE_6_1_1_DIR \
		$VERACODE_CSV \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TARGET_BED
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq )
do
	CREATE_SAMPLE_INFO_ARRAY
	SELECT_PASSING_VARIANTS_PER_SAMPLE
	# BGZIP_AND_TABIX_SAMPLE_VCF
	PASSING_VARIANTS_ON_TARGET_BY_SAMPLE
	PASSING_SNVS_ON_BAIT_BY_SAMPLE
	PASSING_SNVS_ON_TARGET_BY_SAMPLE
	PASSING_INDELS_ON_BAIT_BY_SAMPLE
	PASSING_INDELS_ON_TARGET_BY_SAMPLE
	PASSING_SNVS_TITV_ALL
	TITV_ALL
	PASSING_SNVS_TITV_KNOWN
	TITV_KNOWN
	PASSING_SNVS_TITV_NOVEL
	TITV_NOVEL
	CONCORDANCE_ON_TARGET_PER_SAMPLE
done

##########################################################################
##### BREAKOUTS FOR VARIANT SUMMARY STATS ################################
##########################################################################

GENERATE_STUDY_HAPMAP_SAMPLE_LISTS () 
{
	HAP_MAP_SAMPLE_LIST=(`echo $CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/'$PREFIX'_hapmap_samples.list'`)
	
	MENDEL_SAMPLE_LIST=(`echo $CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/'$PREFIX'_study_samples.list'`)
	
	echo \
		qsub \
			-S /bin/bash \
 			-cwd \
 			-V \
 			-q $QUEUE_LIST \
 			-p $PRIORITY \
 			-j y \
		-N J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
			-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS.log' \
 			-hold_jid I09_VARIANT_ANNOTATOR_REFINED_$PROJECT_MS \
		$SCRIPT_DIR/J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS.sh \
			$CORE_PATH \
			$PROJECT_MS \
			$PREFIX
}


SELECT_SNVS_ALL () 
{
	echo \
	 qsub \
		-S /bin/bash \
 		-cwd \
 		-V \
 		-q $QUEUE_LIST \
 		-p $PRIORITY \
 		-j y \
	-N J10B_SELECT_SNPS_FOR_ALL_SAMPLES_$PROJECT_MS \
	 	-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10B_SELECT_SNPS_FOR_ALL_SAMPLES.log' \
	 	-hold_jid J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/J10B_SELECT_ALL_SAMPLES_SNP.sh \
	 	$JAVA_1_8 \
	 	$GATK_DIR \
	 	$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX
}

SELECT_PASS_STUDY_ONLY_SNP () 
{
	echo \
	qsub \
		-S /bin/bash \
 		-cwd \
 		-V \
 		-q $QUEUE_LIST \
 		-p $PRIORITY \
 		-j y \
	-N J10C_SELECT_PASS_STUDY_ONLY_SNP_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10C_SELECT_PASS_STUDY_ONLY_SNP.log' \	
		-hold_jid J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/J10C_SELECT_PASS_STUDY_ONLY_SNP.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX \
	 	$HAP_MAP_SAMPLE_LIST
}

SELECT_PASS_HAPMAP_ONLY_SNP ()
{
	echo \
	qsub \
		-S /bin/bash \
 		-cwd \
 		-V \
 		-q $QUEUE_LIST \
 		-p $PRIORITY \
 		-j y \
	-N J10D_SELECT_PASS_HAPMAP_ONLY_SNP_$PROJECT_MS \
	 	-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10D_SELECT_PASS_HAPMAP_ONLY_SNP.log' \	
		-hold_jid J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/J10D_SELECT_PASS_HAPMAP_ONLY_SNP.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX \
	 	$MENDEL_SAMPLE_LIST
}

SELECT_INDELS_ALL ()
{
	echo \
	qsub \
		-S /bin/bash \
 		-cwd \
 		-V \
 		-q $QUEUE_LIST \
 		-p $PRIORITY \
 		-j y \
	-N J10E_SELECT_INDELS_FOR_ALL_SAMPLES_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10E_SELECT_INDELS_FOR_ALL_SAMPLES.log' \	
		-hold_jid J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/J10E_SELECT_ALL_SAMPLES_INDELS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX
}

SELECT_PASS_STUDY_ONLY_INDELS ()
{
	echo \
	qsub \
		-S /bin/bash \
 		-cwd \
 		-V \
 		-q $QUEUE_LIST \
 		-p $PRIORITY \
 		-j y \
	-N J10F_SELECT_PASS_STUDY_ONLY_INDEL_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10F_SELECT_PASS_STUDY_ONLY_INDEL.log' \	
		-hold_jid J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/J10F_SELECT_PASS_STUDY_ONLY_INDEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX \
	 	$HAP_MAP_SAMPLE_LIST
}

SELECT_PASS_HAPMAP_ONLY_INDELS ()
{
	echo \
	qsub \
		-S /bin/bash \
 		-cwd \
 		-V \
 		-q $QUEUE_LIST \
 		-p $PRIORITY \
 		-j y \
	-N J10G_SELECT_PASS_HAPMAP_ONLY_INDEL_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10G_SELECT_PASS_HAPMAP_ONLY_INDEL.log' \
		-hold_jid J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/J10G_SELECT_PASS_HAPMAP_ONLY_INDEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX \
	 	$MENDEL_SAMPLE_LIST
}


SELECT_SNVS_ALL_PASS () 
{
	echo \
	qsub \
		-S /bin/bash \
 		-cwd \
 		-V \
 		-q $QUEUE_LIST \
 		-p $PRIORITY \
 		-j y \
	-N J10H_SELECT_SNP_FOR_ALL_SAMPLES_PASS_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10H_SELECT_SNP_FOR_ALL_SAMPLES_PASS.log' \	
		-hold_jid J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/J10H_SELECT_ALL_SAMPLES_SNP_PASS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

SELECT_INDEL_ALL_PASS () 
{
	echo \
	qsub \
		-S /bin/bash \
 		-cwd \
 		-V \
 		-q $QUEUE_LIST \
 		-p $PRIORITY \
 		-j y \
	-N J10I_SELECT_INDEL_FOR_ALL_SAMPLES_PASS_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_J10H_SELECT_INDEL_FOR_ALL_SAMPLES_PASS.log' \		
		-hold_jid J10_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/J10I_SELECT_ALL_SAMPLES_INDEL_PASS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

GENERATE_STUDY_HAPMAP_SAMPLE_LISTS
SELECT_SNVS_ALL
SELECT_PASS_STUDY_ONLY_SNP
SELECT_PASS_HAPMAP_ONLY_SNP
SELECT_INDELS_ALL
SELECT_PASS_STUDY_ONLY_INDELS
SELECT_PASS_HAPMAP_ONLY_INDELS
SELECT_SNVS_ALL_PASS
SELECT_INDEL_ALL_PASS

##########################################################################
######################End of Functions####################################
##########################################################################
