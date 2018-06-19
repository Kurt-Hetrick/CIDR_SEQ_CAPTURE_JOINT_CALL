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
CIDRSEQSUITE_JAVA_DIR="/mnt/linuxtools/JAVA/jre1.7.0_45/bin"
CIDRSEQSUITE_6_1_1_DIR="/mnt/linuxtools/CIDRSEQSUITE/6.1.1"
CIDRSEQSUITE_ANNOVAR_JAVA="/mnt/linuxtools/JAVA/jre1.6.0_25"
CIDRSEQSUITE_DIR_4_0="/mnt/research/tools/LINUX/CIDRSEQSUITE/Version_4_0"
# cp -p /u01/home/hling/cidrseqsuite.props.HGMD /mnt/research/tools/LINUX/00_GIT_REPO_KURT/CIDR_SEQ_CAPTURE_JOINT_CALL/STD_VQSR/cidrseqsuite.props
# 14 June 2018
CIDRSEQSUITE_PROPS_DIR="/mnt/research/tools/LINUX/00_GIT_REPO_KURT/CIDR_SEQ_CAPTURE_JOINT_CALL/STD_VQSR"
CIDRSEQSUITE_7_5_0_DIR="/mnt/research/tools/LINUX/CIDRSEQSUITE/7.5.0"
LAB_QC_DIR="/mnt/linuxtools/CUSTOM_CIDR/EnhancedSequencingQCReport/0.0.2"
	# Copied from \\isilon-cifs\sequencing\CIDRSeqSuiteSoftware\RELEASES\7.0.0\QC_REPORT\EnhancedSequencingQCReport.jar

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

##################################################
##################################################
##### JOINT CALLING PROJECT SET-UP ###############
### WHERE THE MULTI-SAMPLE VCF GETS WRITTEN TO ###
##################################################
##################################################

## This checks to see if bed file directory has been created from a previous run.
## If so, remove it to not interfere with current run

if [ -d $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT ]
then
	rm -rf $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT
fi

############################################################################################
##### MAKE THE FOLLOWING FOLDERS IN THE PROJECT WHERE THE MULTI-SAMPLE VCF IS GOING TO #####
############################################################################################

mkdir -p $CORE_PATH/$PROJECT_MS/{LOGS,COMMAND_LINES}
mkdir -p $CORE_PATH/$PROJECT_MS/TEMP/{BED_FILE_SPLIT,AGGREGATE}
mkdir -p $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/
mkdir -p $CORE_PATH/$PROJECT_MS/GVCF/AGGREGATE
mkdir -p $CORE_PATH/$PROJECT_MS/REPORTS/{ANNOVAR,LAB_PREP_REPORTS_MS}
mkdir -p $CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX

##################################################
### FUNCTIONS FOR JOINT CALLING PROJECT SET-UP ###
##################################################

# grab the reference genome file, dbsnp file and bait bed file for the "project"
## !should do a check here to make sure that there is only one record...!

CREATE_PROJECT_INFO_ARRAY ()
{
PROJECT_INFO_ARRAY=(`sed 's/\r//g' $SAMPLE_SHEET \
	| awk 'BEGIN{FS=","} NR>1 {print $12,$18,$16}' \
	| sed 's/,/\t/g' \
	| sort -k 1,1 \
	| awk '{print $1,$2,$3}' \
	| sort \
	| uniq`)

REF_GENOME=${PROJECT_INFO_ARRAY[0]} # field 12 from the sample sheet
PROJECT_DBSNP=${PROJECT_INFO_ARRAY[1]} # field 18 from the sample sheet
PROJECT_BAIT_BED=${PROJECT_INFO_ARRAY[2]} # field 16 from the sample sheet
}

# GET RID OF ALL THE COMMON BED FILE EFF-UPS,

FORMAT_AND_SCATTER_BAIT_BED() 
{
BED_FILE_PREFIX=(`echo SPLITTED_BED_FILE_`)

# make sure that there is EOF
# remove CARRIAGE RETURNS
# remove CHR PREFIXES (THIS IS FOR GRCH37)
# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
awk 1 $PROJECT_BAIT_BED | sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
>| $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed

# SORT TO GRCH37 ORDER
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

split -l $INTERVALS_DIVIDED -a 4 -d \
$CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_AND_SORTED_BED_FILE.bed \
$CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/$BED_FILE_PREFIX

# ADD A .bed suffix to all of the now splitted files
ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/$BED_FILE_PREFIX* | awk '{print "mv",$0,$0".bed"}' | bash
}

# take all of the project/sample combos in the sample sheet and write a g.vcf file path to a *list file
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

# Run Ben's EnhancedSequencingQCReport which; 
# Generates a QC report for lab specific metrics including Physique Report, Samples Table, Sequencer XML data, Pca and Phoenix.
# Does not check if samples are dropped. 

RUN_LAB_PREP_METRICS ()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N A02-LAB_PREP_METRICS"_"$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PROJECT_MS"-LAB_PREP_METRICS.log" \
	$SCRIPT_DIR/A02_LAB_PREP_METRICS.sh \
		$JAVA_1_8 \
		$LAB_QC_DIR \
		$CORE_PATH \
		$PROJECT_MS \
		$SAMPLE_SHEET
}

############################################################
##### CALL THE ABOVE FUNCTIONS TO SET-UP JOINT CALLING #####
############################################################

CREATE_PROJECT_INFO_ARRAY
FORMAT_AND_SCATTER_BAIT_BED
CREATE_GVCF_LIST
RUN_LAB_PREP_METRICS

# need to add something that will generate lab prep qc metrics

#######################################################################
#######################################################################
################# Scatter of Joint Calling ############################
#######################################################################
#######################################################################

# aggregate all of individual g.vcf into one cohort g.vcf per bed file chunk

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

# genotype the cohort g.vcf chunks

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
	-N B01_GENOTYPE_GVCF_$PROJECT_MS'_'$BED_FILE_NAME \
		-o $CORE_PATH/$PROJECT_MS/LOGS/B01_GENOTYPE_GVCF_$BED_FILE_NAME.log \
		-hold_jid A01_COMBINE_GVCF_$PROJECT_MS'_'$BED_FILE_NAME \
	$SCRIPT_DIR/B01_GENOTYPE_GVCF.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX \
		$BED_FILE_NAME
}

# add dnsnp ID, genotype summaries, gc percentage, variant class, tandem repeat units and homopolymer runs to genotyped g.vcf chunks.

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
	-N C01_VARIANT_ANNOTATOR_$PROJECT_MS'_'$BED_FILE_NAME \
		-o $CORE_PATH/$PROJECT_MS/LOGS/C01_VARIANT_ANNOTATOR_$BED_FILE_NAME.log \
		-hold_jid B01_GENOTYPE_GVCF_$PROJECT_MS'_'$BED_FILE_NAME \
	$SCRIPT_DIR/C01_VARIANT_ANNOTATOR.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX \
		$BED_FILE_NAME \
		$PROJECT_DBSNP
}

# build a string of job names (comma delim) from the variant annotator scatter to store as variable to use as
# hold_jid for the cat variants gather (it's used in the next section after the for loop below)

GENERATE_CAT_VARIANTS_HOLD_ID()
{
	CAT_VARIANTS_HOLD_ID=$CAT_VARIANTS_HOLD_ID'C01_VARIANT_ANNOTATOR_'$PROJECT_MS'_'$BED_FILE_NAME','
}

# for each chunk of the original bed file, do combine gvcfs, then genotype gvcfs, then variant annotator
# then generate a string of all the variant annotator job names submitted

for BED_FILE in $(ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/SPLITTED_BED_FILE*);
do
	BED_FILE_NAME=$(basename $BED_FILE .bed)
	COMBINE_GVCF
	GENOTYPE_GVCF
	VARIANT_ANNOTATOR
	GENERATE_CAT_VARIANTS_HOLD_ID
done

#########################################################
#########################################################
##### VCF Gather and  Genotype Refinement Functions #####
#########################################################
#########################################################

# use cat variants to gather up all of the vcf files above into one big file
# MIGHT WANT TO LOOK INTO GatherVcfs (Picard) here
# Other possibility is MergeVcfs (Picard)...GatherVcfs is supposedly used for scatter operations so hopefully more efficient
# The way that CatVariants is constructed, I think would cause a upper limit to the scatter operation.

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
	-N D01_CAT_VARIANTS_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/D01_CAT_VARIANTS.log \
		-hold_jid $CAT_VARIANTS_HOLD_ID \
	$SCRIPT_DIR/D01_CAT_VARIANTS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

# run the snp vqsr model
# to do: find a better to push out an R version to build the plots
# right now, it's buried inside the shell script itself {grrrr}

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
	-N E01_VARIANT_RECALIBRATOR_SNV_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/E01_VARIANT_RECALIBRATOR_SNV.log \
		-hold_jid D01_CAT_VARIANTS_$PROJECT_MS \
	$SCRIPT_DIR/E01_VARIANT_RECALIBRATOR_SNV.sh \
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

# run the indel vqsr model (concurrently done with the snp model above)
# to do: find a better to push out an R version to build the plots
# right now, it's buried inside the shell script itself {grrrr}

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
	-N E02_VARIANT_RECALIBRATOR_INDEL_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/E02_VARIANT_RECALIBRATOR_INDEL.log \
		-hold_jid D01_CAT_VARIANTS_$PROJECT_MS \
	$SCRIPT_DIR/E02_VARIANT_RECALIBRATOR_INDEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$ONEKG_INDELS_VCF \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

# apply the snp vqsr model to the full vcf
# this wait for both the snp and indel models to be done generating before running.

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
	-N F01_APPLY_RECALIBRATION_SNV_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/F01_APPLY_RECALIBRATION_SNV.log \
		-hold_jid E01_VARIANT_RECALIBRATOR_SNV_$PROJECT_MS \
	$SCRIPT_DIR/F01_APPLY_RECALIBRATION_SNV.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

# now apply the indel vqsr model to the full vcf file

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
	-N G01_APPLY_RECALIBRATION_INDEL_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/G01_APPLY_RECALIBRATION_INDEL.log \
		-hold_jid F01_APPLY_RECALIBRATION_SNV_$PROJECT_MS,\
 E02_VARIANT_RECALIBRATOR_INDEL_$PROJECT_MS \
	$SCRIPT_DIR/G01_APPLY_RECALIBRATION_INDEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

CAT_VARIANTS
VARIANT_RECALIBRATOR_SNV
VARIANT_RECALIBRATOR_INDEL
APPLY_RECALIBRATION_SNV
APPLY_RECALIBRATION_INDEL

##################################################
##################################################
##### SCATTER FOR GENOTYPE REFINEMENT ############
##################################################
##################################################

# do a scatter of genotype refinement using the same chunked bed files use to the g.vcf aggregation
# external priors used are the final 1kg genomes dataset, exac v0.3, no family priors used (no ped file)

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
	-N H01_CALCULATE_GENOTYPE_POSTERIORS_$PROJECT_MS"_"$BED_FILE_NAME \
		-o $CORE_PATH/$PROJECT_MS/LOGS/H01_CALCULATE_GENOTYPE_POSTERIORS_$BED_FILE_NAME".log" \
		-hold_jid G01_APPLY_RECALIBRATION_INDEL_$PROJECT_MS \
	$SCRIPT_DIR/H01_CALCULATE_GENOTYPE_POSTERIORS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$P3_1KG \
		$ExAC \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX \
		$BED_FILE_NAME
}

# recalculate the genotype summaries for the now refined genotypes for each vcf chunk

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
	-N I01_VARIANT_ANNOTATOR_REFINED_$PROJECT_MS"_"$BED_FILE_NAME \
		-o $CORE_PATH/$PROJECT_MS/LOGS/I01_VARIANT_ANNOTATOR_REFINED_$BED_FILE_NAME".log" \
		-hold_jid H01_CALCULATE_GENOTYPE_POSTERIORS_$PROJECT_MS"_"$BED_FILE_NAME \
	$SCRIPT_DIR/I01_VARIANT_ANNOTATOR_REFINED.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$PROJECT_DBSNP \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX \
		$BED_FILE_NAME
}

GENERATE_CAT_REFINED_VARIANTS_HOLD_ID()
{
	CAT_REFINED_VARIANTS_HOLD_ID=$CAT_REFINED_VARIANTS_HOLD_ID'I01_VARIANT_ANNOTATOR_REFINED_'$PROJECT_MS'_'$BED_FILE_NAME','
}

# for each chunk of the original bed file, do combine gvcfs, then genotype gvcfs, then variant annotator
# then generate a string of all the variant annotator job names submitted

for BED_FILE in $(ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/SPLITTED_BED_FILE*);
do
	BED_FILE_NAME=$(basename $BED_FILE .bed)
	CALCULATE_GENOTYPE_POSTERIORS
	VARIANT_ANNOTATOR_REFINED
	GENERATE_CAT_VARIANTS_HOLD_ID
done

#########################################################
#########################################################
##### GT Refined VCF Gather #############################
##### Multi-Sample VCF ANNOVAR ##########################
##### VARIANT SUMMARY STATS VCF BREAKOUTS ###############
#########################################################
#########################################################

# use cat variants to gather up all of the gt refined, reannotated vcf files above into one big file

CAT_REFINED_VARIANTS()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N J01_CAT_REFINED_VARIANTS_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/J01_CAT_REFINED_VARIANTS.log \
		-hold_jid $CAT_REFINED_VARIANTS_HOLD_ID \
	$SCRIPT_DIR/J01_CAT_REFINED_VARIANTS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

# run annovar on the final gt refined vcf file
RUN_ANNOVAR()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST",bigmem.q" \
		-p $PRIORITY \
		-j y \
		-pe slots 5 \
		-R y \
		-l mem_free=300G \
	-N K01_ANNOVAR_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/K01_ANNOVAR.log \
		-hold_jid J01_CAT_REFINED_VARIANTS_$PROJECT_MS \
	$SCRIPT_DIR/K01_ANNOVAR.sh \
		$CIDRSEQSUITE_ANNOVAR_JAVA \
		$CIDRSEQSUITE_DIR_4_0 \
		$CIDRSEQSUITE_PROPS_DIR \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

#################################################################################################
### generate separate sample lists for hapmap samples and study samples #########################
### these are to do breakouts of the refined multi-sample vcf for Hua's variant summary stats ###
#################################################################################################

# generate list files by parsing the header of the final ms vcf file
GENERATE_STUDY_HAPMAP_SAMPLE_LISTS () 
{
	HAP_MAP_SAMPLE_LIST=(`echo $CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/'$PREFIX'_hapmap_samples.list'`)
	
	MENDEL_SAMPLE_LIST=(`echo $CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/'$PREFIX'_study_samples.list'`)
	
	# technically don't have to wait on the gather to happen to do this, but for simplicity sake...
	# if performance becomes an issue then can revisit
	
	echo \
		qsub \
			-S /bin/bash \
 			-cwd \
 			-V \
 			-q $QUEUE_LIST \
 			-p $PRIORITY \
 			-j y \
		-N K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
			-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS.log' \
 			-hold_jid J01_CAT_REFINED_VARIANTS_$PROJECT_MS \
		$SCRIPT_DIR/K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS.sh \
			$CORE_PATH \
			$PROJECT_MS \
			$PREFIX
}

# select all the snp sites
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
	-N K02A01_SELECT_SNPS_FOR_ALL_SAMPLES_$PROJECT_MS \
	 	-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02A01_SELECT_SNPS_FOR_ALL_SAMPLES.log' \
	 	-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/K02A01_SELECT_ALL_SAMPLES_SNP.sh \
	 	$JAVA_1_8 \
	 	$GATK_DIR \
	 	$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX
}

# select only passing snp sites that are polymorphic for the study samples
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
	-N K02A02_SELECT_PASS_STUDY_ONLY_SNP_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02A02_SELECT_PASS_STUDY_ONLY_SNP.log' \	
		-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/K02A02_SELECT_PASS_STUDY_ONLY_SNP.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX \
	 	$HAP_MAP_SAMPLE_LIST
}

# select only passing snp sites that are polymorphic for the hapmap samples
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
	-N K02A03_SELECT_PASS_HAPMAP_ONLY_SNP_$PROJECT_MS \
	 	-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02A03_SELECT_PASS_HAPMAP_ONLY_SNP.log' \	
		-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/K02A03_SELECT_PASS_HAPMAP_ONLY_SNP.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX \
	 	$MENDEL_SAMPLE_LIST
}

# select all the indel (and mixed) sites
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
	-N K02A04_SELECT_INDELS_FOR_ALL_SAMPLES_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02A04_SELECT_INDELS_FOR_ALL_SAMPLES.log' \	
		-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/K02A04_SELECT_ALL_SAMPLES_INDELS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX
}

# select only passing indel/mixed sites that are polymorphic for the study samples
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
	-N K02A05_SELECT_PASS_STUDY_ONLY_INDEL_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02A05_SELECT_PASS_STUDY_ONLY_INDEL.log' \	
		-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/K02A05_SELECT_PASS_STUDY_ONLY_INDEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX \
	 	$HAP_MAP_SAMPLE_LIST
}

# select only passing indel/mixed sites that are polymorphic for the hapmap samples
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
	-N K02A06_SELECT_PASS_HAPMAP_ONLY_INDEL_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02A06_SELECT_PASS_HAPMAP_ONLY_INDEL.log' \
		-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/K02A06_SELECT_PASS_HAPMAP_ONLY_INDEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
	 	$CORE_PATH \
	 	$PROJECT_MS \
	 	$PREFIX \
	 	$MENDEL_SAMPLE_LIST
}

# select all passing snp sites
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
	-N K02A07_SELECT_SNP_FOR_ALL_SAMPLES_PASS_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02A07_SELECT_SNP_FOR_ALL_SAMPLES_PASS.log' \	
		-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/K02A07_SELECT_ALL_SAMPLES_SNP_PASS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

# select all passing indel/mixed sites
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
	-N K02A08_SELECT_INDEL_FOR_ALL_SAMPLES_PASS_$PROJECT_MS \
		-o $CORE_PATH/$PROJECT_MS/LOGS/$PREFIX'_K02A08_SELECT_INDEL_FOR_ALL_SAMPLES_PASS.log' \		
		-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
	$SCRIPT_DIR/K02A08_SELECT_ALL_SAMPLES_INDEL_PASS.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$REF_GENOME \
		$CORE_PATH \
		$PROJECT_MS \
		$PREFIX
}

CAT_REFINED_VARIANTS
RUN_ANNOVAR
GENERATE_STUDY_HAPMAP_SAMPLE_LISTS
SELECT_SNVS_ALL
SELECT_PASS_STUDY_ONLY_SNP
SELECT_PASS_HAPMAP_ONLY_SNP
SELECT_INDELS_ALL
SELECT_PASS_STUDY_ONLY_INDELS
SELECT_PASS_HAPMAP_ONLY_INDELS
SELECT_SNVS_ALL_PASS
SELECT_INDEL_ALL_PASS

#######################################################################
#######################################################################
################### Start of Sample Breakouts #########################
#######################################################################
#######################################################################

# NEED TO ADD MIXED. REMOVED MIXED FROM INDELS.

# for each unique sample id in the sample sheet grab the bed files, ref genome, project and store as an array
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

# for each sample make a bunch directories if not already present in the samples defined project directory
MAKE_PROJ_DIR_TREE ()
{
	mkdir -p \
	$CORE_PATH/$PROJECT_SAMPLE/{TEMP,LOGS,COMMAND_LINES} \
	$CORE_PATH/$PROJECT_SAMPLE/INDEL/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
	$CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
	$CORE_PATH/$PROJECT_SAMPLE/MIXED/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
	$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
	$CORE_PATH/$PROJECT_SAMPLE/REPORTS/{ANNOVAR,QC_REPORT_PREP_MS,QC_REPORTS,LAB_PREP_REPORTS,TI_TV_MS,CONCORDANCE_MS} \
	$CORE_PATH/$PROJECT_SAMPLE/TEMP/$SM_TAG"_MS_CONCORDANCE"
}

# for each sample, make a vcf containing all passing variants
SELECT_PASSING_VARIANTS_PER_SAMPLE ()
{
	echo \
	qsub \
		-S /bin/bash \
		-cwd \
		-V \
		-q $QUEUE_LIST \
		-p $PRIORITY \
		-j y \
	-N K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_MS/LOGS/K03_SELECT_VARIANTS_FOR_SAMPLE_$SM_TAG.log \
		-hold_jid J01_CAT_REFINED_VARIANTS_$PROJECT_MS \
	$SCRIPT_DIR/K03_SELECT_VARIANTS_FOR_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$PREFIX
}

# for each sample, make a vcf containing all passing variants and only falls in the on target bed file
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
	-N K03A01_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-hold_jid K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A01_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
	$SCRIPT_DIR/K03A01_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TARGET_BED
}

# for each sample use the passing on target snvs to calculate concordance and het sensitivity to array genotypes.
# reconfigure using the new concordance tool.
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
		-N K03A01-1_CONCORDANCE_ON_TARGET_PER_SAMPLE_$UNIQUE_ID_SM_TAG \
			-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A01-1_CONCORDANCE_ON_TARGET_PER_SAMPLE_$SAMPLE.log \
			-hold_jid K03A01_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		$SCRIPT_DIR/K03A01-1_CONCORDANCE_ON_TARGET_PER_SAMPLE.sh \
			$JAVA_1_8 \
			$CIDRSEQSUITE_7_5_0_DIR \
			$VERACODE_CSV \
			$CORE_PATH \
			$PROJECT_SAMPLE \
			$SM_TAG \
			$TARGET_BED
	}

# for each sample, make a vcf containing only passing snvs
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
	-N K03A02_PASSING_SNVS_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A02_PASSING_SNVS_ON_BAIT_BY_SAMPLE_$SAMPLE.log \
		-hold_jid K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/K03A02_PASSING_SNVS_ON_BAIT_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG
}

# for each sample, make a vcf containing only passing snvs that fall within the on target bed file
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
	-N K03A03_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A03_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
		-hold_jid K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/K03A03_PASSING_SNVS_ON_TARGET_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TARGET_BED
}

# for each sample, make a vcf containing only passing indels
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
	-N K03A04_PASSING_INDELS_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A04_PASSING_INDELS_ON_BAIT_BY_SAMPLE_$SAMPLE.log \
		-hold_jid K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/K03A04_PASSING_INDELS_ON_BAIT_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG
}

# for each sample, make a vcf containing only passing indels that fall within the on target bed file
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
	-N K03A05_PASSING_INDELS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A05_PASSING_INDELS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
		-hold_jid K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/K03A05_PASSING_INDELS_ON_TARGET_BY_SAMPLE.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TARGET_BED
}

# for each sample, make a vcf containing only passing snvs that fall within the on titv bed file
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
	-N K03A06_PASSING_SNVS_TITV_ALL_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A06_PASSING_SNVS_TITV_ALL_$SAMPLE.log \
		-hold_jid K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/K03A06_PASSING_SNVS_TITV_ALL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TITV_BED
}

	# for each sample, calculate ti/tv for all passing snvs within the ti/tv bed file
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
		-N K03A06-1_TITV_ALL_$UNIQUE_ID_SM_TAG \
			-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A06-1_TITV_ALL_$SAMPLE.log \
			-hold_jid K03_PASSING_SNVS_TITV_ALL_$UNIQUE_ID_SM_TAG \
		$SCRIPT_DIR/K03A06-1_TITV_ALL.sh \
			$SAMTOOLS_0118_DIR \
			$CORE_PATH \
			$PROJECT_SAMPLE \
			$SM_TAG
	}

# for each sample, make a vcf containing only passing snvs that fall within the on titv bed file and are in dbsnp129
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
	-N K03A07_PASSING_SNVS_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A07_PASSING_SNVS_TITV_KNOWN_$SAMPLE.log \
		-hold_jid K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/K03A07_PASSING_SNVS_TITV_KNOWN.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$KNOWN_SNPS \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TITV_BED
}

	# for each sample, calculate ti/tv for all passing snvs within the ti/tv bed file that are in dbsnp129
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
		-N K03A07-1_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
			-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A07-1_TITV_KNOWN_$SAMPLE.log \
			-hold_jid K03A07_PASSING_SNVS_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
		$SCRIPT_DIR/K03A07-1_TITV_KNOWN.sh \
			$SAMTOOLS_0118_DIR \
			$CORE_PATH \
			$PROJECT_SAMPLE \
			$SM_TAG
	}

# for each sample, make a vcf containing only passing snvs that fall within the on titv bed file and are NOT in dbsnp129
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
	-N K03A08_PASSING_SNVS_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
		-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A08_PASSING_SNVS_TITV_NOVEL_$SAMPLE.log \
		-hold_jid K03_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
	$SCRIPT_DIR/K03A08_PASSING_SNVS_TITV_NOVEL.sh \
		$JAVA_1_8 \
		$GATK_DIR \
		$SAMPLE_REF_GENOME \
		$KNOWN_SNPS \
		$CORE_PATH \
		$PROJECT_SAMPLE \
		$SM_TAG \
		$TITV_BED
}

	# for each sample, calculate ti/tv for all passing snvs within the ti/tv bed file that are NOT in dbsnp129
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
		-N K03A08-1_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
			-o $CORE_PATH/$PROJECT_SAMPLE/LOGS/K03A08-1_TITV_NOVEL_$SAMPLE.log \
			-hold_jid K03A08_PASSING_SNVS_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
		$SCRIPT_DIR/K03A08-1_TITV_NOVEL.sh \
			$SAMTOOLS_0118_DIR \
			$CORE_PATH \
			$PROJECT_SAMPLE \
			$SM_TAG
	}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq )
do
	CREATE_SAMPLE_INFO_ARRAY
	MAKE_PROJ_DIR_TREE
	SELECT_PASSING_VARIANTS_PER_SAMPLE
	PASSING_VARIANTS_ON_TARGET_BY_SAMPLE
	PASSING_SNVS_ON_BAIT_BY_SAMPLE
	PASSING_SNVS_ON_TARGET_BY_SAMPLE
	CONCORDANCE_ON_TARGET_PER_SAMPLE
	PASSING_INDELS_ON_BAIT_BY_SAMPLE
	PASSING_INDELS_ON_TARGET_BY_SAMPLE
	PASSING_SNVS_TITV_ALL
	TITV_ALL
	PASSING_SNVS_TITV_KNOWN
	TITV_KNOWN
	PASSING_SNVS_TITV_NOVEL
	TITV_NOVEL
done

##########################################################################
######################End of Functions####################################
##########################################################################

# need to create qc reports, aneuploidy reports and per chr verifybamid reports for the release
