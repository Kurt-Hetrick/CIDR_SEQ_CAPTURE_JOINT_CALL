#! /bin/bash

###################
# INPUT VARIABLES #
###################

	PROJECT_MS=$1 # the project where the multi-sample vcf is being written to
	SAMPLE_SHEET=$2 # full/relative path to the sample sheet
	PREFIX=$3 # prefix name that you want to give the multi-sample vcf
	PRIORITY=$4 # default is -13. do not supply this argument unless you want to change from the default. range is -1 to -1023.

		if [[ ! $PRIORITY ]]
			then
			PRIORITY="-13"
		fi

	NUMBER_OF_BED_FILES=$5 # scatter count, if not supplied then the default is what is below. If you want to change this is you have to supply an input for priority as well.

		if [[ ! $NUMBER_OF_BED_FILES ]]
			then
			NUMBER_OF_BED_FILES=500
		fi

###########################
# CORE VARIABLES/SETTINGS #
###########################

	# gcc is so that it can be pushed out to the compute nodes via qsub (-V)
	module load gcc/7.2.0

	# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED
	SCRIPT_DIR="/mnt/research/tools/LINUX/00_GIT_REPO_KURT/CIDR_SEQ_CAPTURE_JOINT_CALL/HARD_FILTER_ALL"

	# Directory where sequencing projects are located
	CORE_PATH="/mnt/research/active"

	# Generate a list of active queue and remove the ones that I don't want to use
	QUEUE_LIST=`qstat -f -s r \
		| egrep -v "^[0-9]|^-|^queue|^ " \
		| cut -d @ -f 1 \
		| sort \
		| uniq \
		| egrep -v "bigmem.q|all.q|cgc.q|programmers.q|rhel7.q|qtest.q|bigdata.q|uhoh.q|testcgc.q" \
		| datamash collapse 1 \
		| awk '{print $1}'`

	 # Use bigmem.q for ANNOVAR in addition to everything else
	ANNOVAR_QUEUE_LIST=`qstat -f -s r \
		| egrep -v "^[0-9]|^-|^queue|^ " \
		| cut -d @ -f 1 \
		| sort \
		| uniq \
		| egrep -v "all.q|cgc.q|programmers.q|qtest.q|uhoh.q" \
		| datamash collapse 1 \
		| awk '{print $1}'`

	 # | awk '{print $1,"-l \x27hostname=!DellR730-03\x27"}'`

	# EVENTUALLY I WANT THIS SET UP AS AN OPTION WITH A DEFAULT OF X

	PRIORITY="-40"

	# eventually, i want to push this out to something...maybe in the vcf file header.
	PIPELINE_VERSION=`git --git-dir=$SCRIPT_DIR/../.git --work-tree=$SCRIPT_DIR/.. log --pretty=format:'%h' -n 1`

	# generate a random number b/w "0 - 32767" to be used for the job name for variant annotator both pre and post refinement
	# this is to help cut down on the job name length so I can increase the scatter count
	HACK=(`echo $RANDOM`)

	# explicitly setting this b/c not everybody has had the $HOME directory transferred and I'm not going to through
	# and figure out who does and does not have this set correctly
	umask 0007

	# TIMESTAMP=`date '+%F.%H-%M-%S'`

#####################
# PIPELINE PROGRAMS #
#####################

	JAVA_1_8="/mnt/linuxtools/JAVA/jdk1.8.0_73/bin"
	BEDTOOLS_DIR="/mnt/linuxtools/BEDTOOLS/bedtools-2.22.0/bin"
	GATK_DIR="/mnt/linuxtools/GATK/GenomeAnalysisTK-3.7"
	SAMTOOLS_0118_DIR="/mnt/linuxtools/SAMTOOLS/samtools-0.1.18"
		# for ti/tv...vcftools.pl vcfstats
	TABIX_DIR="/mnt/linuxtools/TABIX/tabix-0.2.6"
	CIDRSEQSUITE_JAVA_DIR="/mnt/linuxtools/JAVA/jre1.7.0_45/bin"
	CIDRSEQSUITE_6_1_1_DIR="/mnt/linuxtools/CIDRSEQSUITE/6.1.1"
	CIDRSEQSUITE_ANNOVAR_JAVA="/mnt/linuxtools/JAVA/jdk1.8.0_73/bin"
	CIDRSEQSUITE_DIR_4_0="/mnt/research/tools/LINUX/CIDRSEQSUITE/Version_4_0"
	CIDRSEQSUITE_PROPS_DIR="/mnt/research/tools/LINUX/00_GIT_REPO_KURT/CIDR_SEQ_CAPTURE_JOINT_CALL/STD_VQSR"
		# cp -p /u01/home/hling/cidrseqsuite.props.HGMD /mnt/research/tools/LINUX/00_GIT_REPO_KURT/CIDR_SEQ_CAPTURE_JOINT_CALL/STD_VQSR/cidrseqsuite.props
		# 14 June 2018
	CIDRSEQSUITE_7_5_0_DIR="/mnt/research/tools/LINUX/CIDRSEQSUITE/7.5.0"
	LAB_QC_DIR="/mnt/linuxtools/CUSTOM_CIDR/EnhancedSequencingQCReport/0.1.0"
		# Copied from /mnt/research/tools/LINUX/CIDRSEQSUITE/pipeline_dependencies/QC_REPORT/EnhancedSequencingQCReport.jar
		# md5 f979bb4dc8d97113735ef17acd3a766e  EnhancedSequencingQCReport.jar
	SAMTOOLS_DIR="/mnt/linuxtools/ANACONDA/anaconda2-5.0.0.1/bin"
		# This is samtools version 1.7
		# I have no idea why other users other than me cannot index a cram file with a version of samtools that I built from the source
		# Apparently the version that I built with Anaconda works for other users, but it performs REF_CACHE first...
	DATAMASH_DIR="/mnt/linuxtools/DATAMASH/datamash-1.0.6"
	R_DIRECTORY="/mnt/linuxtools/R/R-3.1.1/bin"
	GATK_DIR_4011="/mnt/linuxtools/GATK/gatk-4.0.11.0"
	PICARD_DIR="/mnt/linuxtools/PICARD/picard-2.18.25"

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
	MERGED_CUTTING_BED_FILE="/mnt/research/active/H_Cutting_CFTR_WGHum-SeqCustom_1_Reanalysis/BED_Files/H_Cutting_phase_1plus2_super_file.bed"
	REF_DICT="/mnt/research/tools/PIPELINE_FILES/bwa_mem_0.7.5a_ref/human_g1k_v37_decoy.dict"
	# what is used in the clinical exome pipeline line...at least at first. basically longest transcript per translated gene
	REF_SEQ_TRANSCRIPTS="/mnt/research/tools/PIPELINE_FILES/GRCh37_aux_files/RefSeq.Unique.GRCh37.FINAL.19Feb2018.bed"

##################################################
##################################################
##### JOINT CALLING PROJECT SET-UP ###############
### WHERE THE MULTI-SAMPLE VCF GETS WRITTEN TO ###
##################################################
##################################################

	# This checks to see if bed file directory and split gvcf list has been created from a previous run.
	# If so, remove them to not interfere with current run

	if [ -d $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT ]
		then
			rm -rf $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT
	fi

	if [ -d $CORE_PATH/$PROJECT_MS/TEMP/SPLIT_LIST ]
		then
			rm -rf $CORE_PATH/$PROJECT_MS/TEMP/SPLIT_LIST
	fi

############################################################################################
##### MAKE THE FOLLOWING FOLDERS IN THE PROJECT WHERE THE MULTI-SAMPLE VCF IS GOING TO #####
############################################################################################

	mkdir -p $CORE_PATH/$PROJECT_MS/{LOGS,COMMAND_LINES}
	mkdir -p $CORE_PATH/$PROJECT_MS/TEMP/{BED_FILE_SPLIT,AGGREGATE,QC_REPORT_PREP_$PREFIX,SPLIT_LIST}
	mkdir -p $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/
	mkdir -p $CORE_PATH/$PROJECT_MS/GVCF/AGGREGATE
	mkdir -p $CORE_PATH/$PROJECT_MS/REPORTS/{ANNOVAR,LAB_PREP_REPORTS_MS,QC_REPORTS,QC_REPORT_PREP_$PREFIX}
	mkdir -p $CORE_PATH/$PROJECT_MS/TEMP/ANNOVAR/$PREFIX
	mkdir -p $CORE_PATH/$PROJECT_MS/LOGS/{A01_COMBINE_GVCF,B01_GENOTYPE_GVCF,C01_VARIANT_ANNOTATOR,H01_CALCULATE_GENOTYPE_POSTERIORS,I01_VARIANT_ANNOTATOR_REFINED}

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

				if [[ $PROJECT_MS = "H_Cutting"* ]];
					then
						PROJECT_BAIT_BED=${MERGED_CUTTING_BED_FILE}
					else
						PROJECT_BAIT_BED=${PROJECT_BAIT_BED}
				fi

		}

		# Keep this in here to reference...I'll probably going to use this somewhere.
		# generates a chrosome list from a bait bed file at the project level. so I only have to make one iteration per project instead of one per sample

			# CREATE_CHROMOSOME_LIST ()
			# {
			# 	REF_CHROM_SET=$(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $MERGED_MENDEL_BED_FILE \
			# 		| sed -r 's/[[:space:]]+/\t/g' \
			# 		| cut -f 1 \
			# 		| sort \
			# 		| uniq \
			# 		| $DATAMASH_DIR/datamash collapse 1 \
			# 		| sed 's/,/ /g')
			# }

	# get the full path of the last gvcf list file.
	# take the sample sheet and create a gvcf list from that
	# append the two outputs and write PROJECT_MS-#line_count".samples.list"
	# split them into 300 sample chunks

		CREATE_GVCF_LIST ()
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
			 | egrep -v "idx|tbi|md5" \
			>| $CORE_PATH'/'$PROJECT_MS'/'$PROJECT_MS'-'$TOTAL_SAMPLES'.samples.gvcf.list'

			GVCF_LIST=(`echo $CORE_PATH'/'$PROJECT_MS'/'$PROJECT_MS'-'$TOTAL_SAMPLES'.samples.gvcf.list'`)

			# Take the list above and split it into groups of 300
				split -l 300 -a 4 -d $GVCF_LIST \
				$CORE_PATH/$PROJECT_MS/TEMP/SPLIT_LIST/

			# append *list suffix to output
				ls $CORE_PATH/$PROJECT_MS/TEMP/SPLIT_LIST/* \
					| awk '{print "mv",$0,$0".list"}' \
					| bash	
		}

	# GET RID OF ALL THE COMMON BED FILE EFF-UPS and break up into chunks,

		FORMAT_AND_SCATTER_BAIT_BED ()
		{
				BED_FILE_PREFIX=(`echo BF`)

				# make sure that there is EOF
				# remove CARRIAGE RETURNS
				# remove CHR PREFIXES (THIS IS FOR GRCH37)
				# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.

					awk 1 $PROJECT_BAIT_BED \
						| sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
					>| $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed

				# SORT TO REF ORDER

					(awk '$1~/^[0-9]/' $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed | sort -k1,1n -k2,2n ; \
						awk '$1=="X"' $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed | sort -k 2,2n ; \
						awk '$1=="Y"' $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_BED_FILE.bed | sort -k 2,2n) \
					>| $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_AND_SORTED_BED_FILE.bed

				# get a line count for the number of for the bed file above
				# divide the line count by the number of mini-bed files you want
				# if there is a remainder round up the next integer

				# this somehow sort of works, but the math ends up being off...
				# if I wanted 1000 fold scatter gather, this would actually create 997.
				# and I don't see how this actually rounds up, but it must otherwise split would not work.

					INTERVALS_DIVIDED=`wc -l $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_AND_SORTED_BED_FILE.bed \
						| awk '{print $1"/""'$NUMBER_OF_BED_FILES'"}' \
						| bc \
						| awk '{print $0+1}'`

					split -l $INTERVALS_DIVIDED \
						-a 4 \
						-d \
					$CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/FORMATTED_AND_SORTED_BED_FILE.bed \
					$CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/$BED_FILE_PREFIX

				# ADD A .bed suffix to all of the now splitted files

					ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/$BED_FILE_PREFIX* \
						| awk '{print "mv",$0,$0".bed"}' \
						| bash
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
	# CREATE_CHROMOSOME_LIST
	FORMAT_AND_SCATTER_BAIT_BED
	CREATE_GVCF_LIST
	RUN_LAB_PREP_METRICS
	echo sleep 0.1s

#######################################################################
#######################################################################
################# Scatter of Joint Calling ############################
#######################################################################
#######################################################################

	# aggregate all of individual g.vcf into one cohort g.vcf per bed file chunk

		COMBINE_GVCF ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
				-j y \
			-N 'A01_COMBINE_GVCF_'$PROJECT_MS'_'$PGVCF_LIST_NAME'_'$BED_FILE_NAME \
				-o $CORE_PATH/$PROJECT_MS/LOGS/A01_COMBINE_GVCF/"A01_COMBINE_GVCF_"$PGVCF_LIST_NAME"_"$BED_FILE_NAME".log" \
			$SCRIPT_DIR/A01_COMBINE_GVCF.sh \
				$JAVA_1_8 \
				$GATK_DIR \
				$REF_GENOME \
				$CORE_PATH \
				$PROJECT_MS \
				$PGVCF_LIST \
				$PREFIX \
				$BED_FILE_NAME
		}

for BED_FILE in $(ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/BF*bed);
	do
		BED_FILE_NAME=$(basename $BED_FILE .bed)
			for PGVCF_LIST in $(ls $CORE_PATH/$PROJECT_MS/TEMP/SPLIT_LIST/*list)
				do
					PGVCF_LIST_NAME=$(basename $PGVCF_LIST .list)
					COMBINE_GVCF
					echo sleep 0.1s
			done
done

#####################################################################

	BUILD_HOLD_ID_GENOTYPE_GVCF ()
	{
		for PROJECT_A in $PROJECT_MS;
		# yeah, so uh, this looks bad, but I just needed a way to set a new project variable that equals the multi-sample project variable.
			do
				GENOTYPE_GVCF_HOLD_ID="-hold_jid "

					for PGVCF_LIST in $(ls $CORE_PATH/$PROJECT_A/TEMP/SPLIT_LIST/*list)
						do
							PGVCF_LIST_NAME=$(basename $PGVCF_LIST .list)
							GENOTYPE_GVCF_HOLD_ID=$GENOTYPE_GVCF_HOLD_ID'A01_COMBINE_GVCF_'$PROJECT_A'_'$PGVCF_LIST_NAME'_'$BED_FILE_NAME','
					done
		done
	}

	# genotype the cohort g.vcf chunks

		GENOTYPE_GVCF ()
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
				-o $CORE_PATH/$PROJECT_MS/LOGS/B01_GENOTYPE_GVCF/B01_GENOTYPE_GVCF_$BED_FILE_NAME.log \
				$GENOTYPE_GVCF_HOLD_ID \
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

		VARIANT_ANNOTATOR ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
				-j y \
			-N C$HACK'_'$BED_FILE_NAME \
				-o $CORE_PATH/$PROJECT_MS/LOGS/C01_VARIANT_ANNOTATOR/C01_VARIANT_ANNOTATOR_$BED_FILE_NAME.log \
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

		GENERATE_CAT_VARIANTS_HOLD_ID ()
		{
			CAT_VARIANTS_HOLD_ID=$CAT_VARIANTS_HOLD_ID'C'$HACK'_'$BED_FILE_NAME','
		}

	# for each chunk of the original bed file, do combine gvcfs, then genotype gvcfs, then variant annotator
	# then generate a string of all the variant annotator job names submitted

for BED_FILE in $(ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/BF*bed);
	do
		BED_FILE_NAME=$(basename $BED_FILE .bed)
		BUILD_HOLD_ID_GENOTYPE_GVCF
		GENOTYPE_GVCF
		echo sleep 0.1s
		VARIANT_ANNOTATOR
		echo sleep 0.1s
		GENERATE_CAT_VARIANTS_HOLD_ID
done

#######################################
#######################################
##### VCF Gather and Hard Filters #####
#######################################
#######################################

	# use cat variants to gather up all of the vcf files above into one big file
	# MIGHT WANT TO LOOK INTO GatherVcfs (Picard) here
	# Other possibility is MergeVcfs (Picard)...GatherVcfs is supposedly used for scatter operations so hopefully more efficient
	# The way that CatVariants is constructed, I think would cause a upper limit to the scatter operation.

		CAT_VARIANTS ()
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

	# Breakout the snps

		EXTRACT_SNV ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N E01_EXTRACT_SNV_$PROJECT_MS \
					-o $CORE_PATH/$PROJECT_MS/LOGS/E01_EXTRACT_SNV.log \
					-hold_jid D01_CAT_VARIANTS_$PROJECT_MS \
				$SCRIPT_DIR/E01_EXTRACT_SNV.sh \
					$JAVA_1_8 \
					$GATK_DIR \
					$REF_GENOME \
					$CORE_PATH \
					$PROJECT_MS \
					$PREFIX
			}

		# APPLY HARD FILTERS TO SNPS

			HARD_FILTER_SNV ()
				{
					echo \
					qsub \
						-S /bin/bash \
						-cwd \
						-V \
						-q $QUEUE_LIST \
						-p $PRIORITY \
						-j y \
					-N E01-A01_FILTER_SNV_$PROJECT_MS \
						-o $CORE_PATH/$PROJECT_MS/LOGS/E01-A01_FILTER_SNV.log \
						-hold_jid E01_EXTRACT_SNV_$PROJECT_MS \
					$SCRIPT_DIR/E01-A01_FILTER_SNV.sh \
						$JAVA_1_8 \
						$GATK_DIR \
						$REF_GENOME \
						$CORE_PATH \
						$PROJECT_MS \
						$PREFIX
				}

	# Breakout the non-snps

		EXTRACT_INDEL_AND_MIXED ()
			{
				echo \
					qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N E02_EXTRACT_INDEL_AND_MIXED_$PROJECT_MS \
					-o $CORE_PATH/$PROJECT_MS/LOGS/E02_EXTRACT_INDEL_AND_MIXED.log \
					-hold_jid D01_CAT_VARIANTS_$PROJECT_MS \
				$SCRIPT_DIR/E02_EXTRACT_INDELS_AND_MIXED.sh \
					$JAVA_1_8 \
					$GATK_DIR \
					$REF_GENOME \
					$CORE_PATH \
					$PROJECT_MS \
					$PREFIX
			}

		# Apply Hard Filters to Indels and Mixed

			HARD_FILTER_INDEL_AND_MIXED ()
				{
					echo \
					qsub \
						-S /bin/bash \
						-cwd \
						-V \
						-q $QUEUE_LIST \
						-p $PRIORITY \
						-j y \
					-N E02-A01_FILTER_INDEL_AND_MIXED_$PROJECT_MS \
						-o $CORE_PATH/$PROJECT_MS/LOGS/E02-A01_FILTER_INDEL_AND_MIXED.log \
						-hold_jid E02_EXTRACT_INDEL_AND_MIXED_$PROJECT_MS \
					$SCRIPT_DIR/E02-A01_FILTER_INDELS_AND_MIXED.sh \
						$JAVA_1_8 \
						$GATK_DIR \
						$REF_GENOME \
						$CORE_PATH \
						$PROJECT_MS \
						$PREFIX
				}

	# COMBINE_VARIANTS

		COMBINE_SNV_INDEL_VARIANTS ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N F01_COMBINE_VARIANTS_$PROJECT_MS \
					-o $CORE_PATH/$PROJECT_MS/LOGS/F01_COMBINE_VARIANTS.log \
					-hold_jid E01-A01_FILTER_SNV_$PROJECT_MS,E02-A01_FILTER_INDEL_AND_MIXED_$PROJECT_MS \
				$SCRIPT_DIR/F01_COMBINE_VARIANTS.sh \
						$JAVA_1_8 \
						$GATK_DIR \
						$REF_GENOME \
						$CORE_PATH \
						$PROJECT_MS \
						$PREFIX
			}


# call cat variants and hard filter variants

	CAT_VARIANTS
	echo sleep 0.1s
	EXTRACT_SNV
	echo sleep 0.1s
	EXTRACT_INDEL_AND_MIXED
	echo sleep 0.1s
	HARD_FILTER_SNV
	echo sleep 0.1s
	HARD_FILTER_INDEL_AND_MIXED
	echo sleep 0.1s
	COMBINE_SNV_INDEL_VARIANTS
	echo sleep 0.1s

##################################################
##################################################
##### SCATTER FOR GENOTYPE REFINEMENT ############
##################################################
##################################################

	# do a scatter of genotype refinement using the same chunked bed files use to the g.vcf aggregation
	# external priors used are the final 1kg genomes dataset, exac v0.3, no family priors used (no ped file)

		CALCULATE_GENOTYPE_POSTERIORS ()
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
				-o $CORE_PATH/$PROJECT_MS/LOGS/H01_CALCULATE_GENOTYPE_POSTERIORS/H01_CALCULATE_GENOTYPE_POSTERIORS_$BED_FILE_NAME".log" \
			-hold_jid F01_COMBINE_VARIANTS_$PROJECT_MS \
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

		VARIANT_ANNOTATOR_REFINED ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
				-j y \
			-N I$HACK"_"$BED_FILE_NAME \
				-o $CORE_PATH/$PROJECT_MS/LOGS/I01_VARIANT_ANNOTATOR_REFINED/I01_VARIANT_ANNOTATOR_REFINED_$BED_FILE_NAME".log" \
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

	# build a string of job names (comma delim) from the variant annotator scatter to store as variable to use as
	# hold_jid for the cat variants gather (it's used in the next section after the for loop below)

		GENERATE_CAT_REFINED_VARIANTS_HOLD_ID ()
		{
			CAT_REFINED_VARIANTS_HOLD_ID=$CAT_REFINED_VARIANTS_HOLD_ID'I'$HACK'_'$BED_FILE_NAME','
		}

	# for each chunk of the original bed file, do calculate_genotype_posteriors, then variant annotator
	# then generate a string of all the variant annotator job names submitted

for BED_FILE in $(ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/BF*);
do
	BED_FILE_NAME=$(basename $BED_FILE .bed)
	CALCULATE_GENOTYPE_POSTERIORS
	echo sleep 0.1s
	VARIANT_ANNOTATOR_REFINED
	echo sleep 0.1s
	GENERATE_CAT_REFINED_VARIANTS_HOLD_ID
done

#########################################################
#########################################################
##### GT Refined VCF Gather #############################
##### Multi-Sample VCF ANNOVAR ##########################
##### VARIANT SUMMARY STATS VCF BREAKOUTS ###############
#########################################################
#########################################################

	# use cat variants to gather up all of the gt refined, reannotated vcf files above into one big file

		CAT_REFINED_VARIANTS ()
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

	# bgzip and index genotype refined vcf

		BGZIP_INDEX_REFINED_VARIANTS ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
				-j y \
			-N J01A01_BGZIP_INDEX_$PROJECT_MS \
				-o $CORE_PATH/$PROJECT_MS/LOGS/J01A01_BGZIP_INDEX_VARIANTS.log \
			-hold_jid J01_CAT_REFINED_VARIANTS_$PROJECT_MS \
			$SCRIPT_DIR/J01A01_BGZIP_INDEX_VARIANTS.sh \
				$TABIX_DIR \
				$CORE_PATH \
				$PROJECT_MS \
				$PREFIX
		}

		# # run annovar on the final gt refined vcf file

			# 	RUN_ANNOVAR ()
			# 	{
			# 		echo \
			# 		qsub \
			# 			-S /bin/bash \
			# 			-cwd \
			# 			-V \
			# 			-q $QUEUE_LIST",bigmem.q" \
			# 			-p $PRIORITY \
			# 			-j y \
			# 			-pe slots 5 \
			# 			-R y \
			# 			-l mem_free=300G \
			# 		-N K01_ANNOVAR_$PROJECT_MS \
			# 			-o $CORE_PATH/$PROJECT_MS/LOGS/K01_ANNOVAR.log \
			# 			-hold_jid J01_CAT_REFINED_VARIANTS_$PROJECT_MS \
			# 		$SCRIPT_DIR/K01_ANNOVAR.sh \
			# 			$CIDRSEQSUITE_ANNOVAR_JAVA \
			# 			$CIDRSEQSUITE_DIR_4_0 \
			# 			$CIDRSEQSUITE_PROPS_DIR \
			# 			$CORE_PATH \
			# 			$PROJECT_MS \
			# 			$PREFIX
			# 	}

# cat refined variants, annovar

	CAT_REFINED_VARIANTS
	echo sleep 0.1s
	# RUN_ANNOVAR
	# echo sleep 0.1s
	BGZIP_INDEX_REFINED_VARIANTS
	echo sleep 0.1s

#########################################################
#########################################################
##### VARIANT SUMMARY STATS VCF BREAKOUTS ###############
#########################################################
#########################################################

	###########################################################################
	### generate separate sample lists for hapmap samples and study samples ###
	### these are to do breakouts for Hua's variant summary stats #############
	###########################################################################

		# generate list files by parsing the header of the final ms vcf file
			GENERATE_STUDY_HAPMAP_SAMPLE_LISTS ()
			{
				HAP_MAP_SAMPLE_LIST=(`echo $CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/'$PREFIX'_hapmap_samples.args'`)

				MENDEL_SAMPLE_LIST=(`echo $CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/'$PREFIX'_study_samples.args'`)

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
						-o $CORE_PATH/$PROJECT_MS/LOGS/'K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS-'$PREFIX'.log' \
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
					-o $CORE_PATH/$PROJECT_MS/LOGS/'K02A01_SELECT_SNPS_FOR_ALL_SAMPLES-'$PREFIX'.log' \
				-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
				$SCRIPT_DIR/K02A01_SELECT_ALL_SAMPLES_SNP.sh \
				 	$JAVA_1_8 \
					$GATK_DIR_4011 \
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
					-o $CORE_PATH/$PROJECT_MS/LOGS/'K02A02_SELECT_PASS_STUDY_ONLY_SNP-'$PREFIX'.log' \
				-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
				$SCRIPT_DIR/K02A02_SELECT_PASS_STUDY_ONLY_SNP.sh \
					$JAVA_1_8 \
					$GATK_DIR_4011 \
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
					-o $CORE_PATH/$PROJECT_MS/LOGS/'K02A03_SELECT_PASS_HAPMAP_ONLY_SNP-'$PREFIX'.log' \
				-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
				$SCRIPT_DIR/K02A03_SELECT_PASS_HAPMAP_ONLY_SNP.sh \
					$JAVA_1_8 \
					$GATK_DIR_4011 \
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
					-o $CORE_PATH/$PROJECT_MS/LOGS/'K02A04_SELECT_INDELS_FOR_ALL_SAMPLES-'$PREFIX'.log' \
				-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
				$SCRIPT_DIR/K02A04_SELECT_ALL_SAMPLES_INDELS.sh \
					$JAVA_1_8 \
					$GATK_DIR_4011 \
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
					-o $CORE_PATH/$PROJECT_MS/LOGS/'K02A05_SELECT_PASS_STUDY_ONLY_INDEL-'$PREFIX'.log' \
				-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
				$SCRIPT_DIR/K02A05_SELECT_PASS_STUDY_ONLY_INDEL.sh \
					$JAVA_1_8 \
					$GATK_DIR_4011 \
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
					-o $CORE_PATH/$PROJECT_MS/LOGS/'K02A06_SELECT_PASS_HAPMAP_ONLY_INDEL-'$PREFIX'.log' \
				-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
				$SCRIPT_DIR/K02A06_SELECT_PASS_HAPMAP_ONLY_INDEL.sh \
					$JAVA_1_8 \
					$GATK_DIR_4011 \
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
					-o $CORE_PATH/$PROJECT_MS/LOGS/'K02A07_SELECT_SNP_FOR_ALL_SAMPLES_PASS-'$PREFIX'.log' \
				-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
				$SCRIPT_DIR/K02A07_SELECT_ALL_SAMPLES_SNP_PASS.sh \
					$JAVA_1_8 \
					$GATK_DIR_4011 \
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
					-o $CORE_PATH/$PROJECT_MS/LOGS/'K02A08_SELECT_INDEL_FOR_ALL_SAMPLES_PASS-'$PREFIX'.log' \
				-hold_jid K02_GENERATE_STUDY_HAPMAP_SAMPLE_LISTS_$PROJECT_MS \
				$SCRIPT_DIR/K02A08_SELECT_ALL_SAMPLES_INDEL_PASS.sh \
					$JAVA_1_8 \
					$GATK_DIR_4011 \
					$REF_GENOME \
					$CORE_PATH \
					$PROJECT_MS \
					$PREFIX
			}

# call variant summary stat vcf breakouts

	GENERATE_STUDY_HAPMAP_SAMPLE_LISTS
	SELECT_SNVS_ALL
	echo sleep 0.1s
	SELECT_PASS_STUDY_ONLY_SNP
	echo sleep 0.1s
	SELECT_PASS_HAPMAP_ONLY_SNP
	echo sleep 0.1s
	SELECT_INDELS_ALL
	echo sleep 0.1s
	SELECT_PASS_STUDY_ONLY_INDELS
	echo sleep 0.1s
	SELECT_PASS_HAPMAP_ONLY_INDELS
	echo sleep 0.1s
	SELECT_SNVS_ALL_PASS
	echo sleep 0.1s
	SELECT_INDEL_ALL_PASS
	echo sleep 0.1s

#######################################################################
#######################################################################
################### Start of Sample Breakouts #########################
#######################################################################
#######################################################################


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
			DBSNP=${SAMPLE_INFO_ARRAY[4]}
			SAMPLE_REF_GENOME=${SAMPLE_INFO_ARRAY[5]}

			UNIQUE_ID_SM_TAG=$(echo $SM_TAG | sed 's/@/_/g') # If there is an @ in the qsub or holdId name it breaks
			BARCODE_2D=$(echo $SM_TAG | awk '{split($1,SM_TAG,/[@-]/); print SM_TAG[2]}') # SM_TAG = RIS_ID[@-]BARCODE_2D
		}

	# for each sample make a bunch directories if not already present in the samples defined project directory
	# shouldn't need this anymore $SM_TAG"_SCATTER"

		MAKE_PROJ_DIR_TREE ()
		{
			mkdir -p \
			$CORE_PATH/$PROJECT_SAMPLE/{TEMP,LOGS,COMMAND_LINES} \
			$CORE_PATH/$PROJECT_SAMPLE/INDEL/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
			$CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
			$CORE_PATH/$PROJECT_SAMPLE/MIXED/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
			$CORE_PATH/$PROJECT_SAMPLE/VCF/RELEASE/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
			$CORE_PATH/$PROJECT_SAMPLE/REPORTS/{TI_TV_MS,CONCORDANCE_MS,ANNOVAR} \
			$CORE_PATH/$PROJECT_SAMPLE/TEMP/{$SM_TAG"_MS_CONCORDANCE",$SM_TAG"_ANNOVAR",$SM_TAG} \
			$CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG
		}


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
			-N K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01_SELECT_VARIANTS_FOR_SAMPLE_$SM_TAG".log" \
			-hold_jid J01_CAT_REFINED_VARIANTS_$PROJECT_MS \
			$SCRIPT_DIR/K01_SELECT_VARIANTS_FOR_SAMPLE.sh \
				$JAVA_1_8 \
				$GATK_DIR_4011 \
				$SAMPLE_REF_GENOME \
				$CORE_PATH \
				$PROJECT_SAMPLE \
				$PROJECT_MS \
				$SM_TAG \
				$PREFIX
		}


		PASSING_VARIANTS_ON_TARGET_BY_SAMPLE ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
				-j y \
			-N K01A01_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
				-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A01_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
			-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
			$SCRIPT_DIR/K01A01_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE.sh \
				$JAVA_1_8 \
				$GATK_DIR \
				$SAMPLE_REF_GENOME \
				$CORE_PATH \
				$PROJECT_SAMPLE \
				$SM_TAG \
				$TARGET_BED
		}

	########################################################################
	### grabbing per sample SNV only vcf files for on bait and on target ###
	########################################################################

			PASSING_SNVS_ON_BAIT_BY_SAMPLE ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A02_PASSING_SNVS_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A02_PASSING_SNVS_ON_BAIT_BY_SAMPLE_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A02_PASSING_SNVS_ON_BAIT_BY_SAMPLE.sh \
					$JAVA_1_8 \
					$GATK_DIR \
					$SAMPLE_REF_GENOME \
					$CORE_PATH \
					$PROJECT_SAMPLE \
					$SM_TAG
			}

		# for each sample, make a vcf containing only passing snvs that fall within the on target bed file

			PASSING_SNVS_ON_TARGET_BY_SAMPLE ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A03_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A03_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A03_PASSING_SNVS_ON_TARGET_BY_SAMPLE.sh \
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
					CONCORDANCE_ON_TARGET_PER_SAMPLE ()
					{
						echo \
						qsub \
							-S /bin/bash \
							-cwd \
							-V \
							-q $QUEUE_LIST \
							-p $PRIORITY \
							-j y \
						-N K01A03-1_CONCORDANCE_ON_TARGET_PER_SAMPLE_$UNIQUE_ID_SM_TAG \
							-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A03-1_CONCORDANCE_ON_TARGET_$SAMPLE.log \
						-hold_jid K01A03_PASSING_SNVS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
						$SCRIPT_DIR/K03A03-1_CONCORDANCE_ON_TARGET_PER_SAMPLE.sh \
							$JAVA_1_8 \
							$CIDRSEQSUITE_7_5_0_DIR \
							$VERACODE_CSV \
							$CORE_PATH \
							$PROJECT_SAMPLE \
							$SM_TAG \
							$TARGET_BED
					}

	##########################################################################
	### grabbing per sample indel only vcf files for on bait and on target ###
	##########################################################################

		# for each sample, make a vcf containing only passing indels

			PASSING_INDELS_ON_BAIT_BY_SAMPLE ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A04_PASSING_INDELS_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A04_PASSING_INDELS_ON_BAIT_BY_SAMPLE_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A04_PASSING_INDELS_ON_BAIT_BY_SAMPLE.sh \
					$JAVA_1_8 \
					$GATK_DIR \
					$SAMPLE_REF_GENOME \
					$CORE_PATH \
					$PROJECT_SAMPLE \
					$SM_TAG
			}

		# for each sample, make a vcf containing only passing indels that fall within the on target bed file

			PASSING_INDELS_ON_TARGET_BY_SAMPLE ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A05_PASSING_INDELS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A05_PASSING_INDELS_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A05_PASSING_INDELS_ON_TARGET_BY_SAMPLE.sh \
					$JAVA_1_8 \
					$GATK_DIR \
					$SAMPLE_REF_GENOME \
					$CORE_PATH \
					$PROJECT_SAMPLE \
					$SM_TAG \
					$TARGET_BED
			}

	########################################################
	### grabbing per sample snv only vcf files for ti/tv ###
	########################################################

		# for each sample, make a vcf containing only passing snvs that fall within the on titv bed file

			PASSING_SNVS_TITV_ALL ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A06_PASSING_SNVS_TITV_ALL_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A06_PASSING_SNVS_TITV_ALL_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A06_PASSING_SNVS_TITV_ALL.sh \
					$JAVA_1_8 \
					$GATK_DIR \
					$SAMPLE_REF_GENOME \
					$CORE_PATH \
					$PROJECT_SAMPLE \
					$SM_TAG \
					$TITV_BED
			}

			# for each sample, calculate ti/tv for all passing snvs within the ti/tv bed file

				TITV_ALL ()
				{
					echo \
					qsub \
						-S /bin/bash \
						-cwd \
						-V \
						-q $QUEUE_LIST \
						-p $PRIORITY \
						-j y \
					-N K01A06-1_TITV_ALL_$UNIQUE_ID_SM_TAG \
						-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A06-1_TITV_ALL_$SAMPLE.log \
					-hold_jid K01A06_PASSING_SNVS_TITV_ALL_$UNIQUE_ID_SM_TAG \
					$SCRIPT_DIR/K01A06-1_TITV_ALL.sh \
						$SAMTOOLS_0118_DIR \
						$CORE_PATH \
						$PROJECT_SAMPLE \
						$SM_TAG
				}

		# for each sample, make a vcf containing only passing snvs that fall within the on titv bed file and are in dbsnp129

			PASSING_SNVS_TITV_KNOWN ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A07_PASSING_SNVS_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A07_PASSING_SNVS_TITV_KNOWN_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A07_PASSING_SNVS_TITV_KNOWN.sh \
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

				TITV_KNOWN ()
				{
					echo \
					qsub \
						-S /bin/bash \
						-cwd \
						-V \
						-q $QUEUE_LIST \
						-p $PRIORITY \
						-j y \
						-N K01A07-1_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
							-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A07-1_TITV_KNOWN_$SAMPLE.log \
						-hold_jid K01A07_PASSING_SNVS_TITV_KNOWN_$UNIQUE_ID_SM_TAG \
						$SCRIPT_DIR/K01A07-1_TITV_KNOWN.sh \
						$SAMTOOLS_0118_DIR \
						$CORE_PATH \
						$PROJECT_SAMPLE \
						$SM_TAG
				}

		# for each sample, make a vcf containing only passing snvs that fall within the on titv bed file and are NOT in dbsnp129

			PASSING_SNVS_TITV_NOVEL ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A08_PASSING_SNVS_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A08_PASSING_SNVS_TITV_NOVEL_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A08_PASSING_SNVS_TITV_NOVEL.sh \
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

				TITV_NOVEL ()
				{
					echo \
					qsub \
						-S /bin/bash \
						-cwd \
						-V \
						-q $QUEUE_LIST \
						-p $PRIORITY \
						-j y \
					-N K01A08-1_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
						-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A08-1_TITV_NOVEL_$SAMPLE.log \
					-hold_jid K01A08_PASSING_SNVS_TITV_NOVEL_$UNIQUE_ID_SM_TAG \
					$SCRIPT_DIR/K01A08-1_TITV_NOVEL.sh \
						$SAMTOOLS_0118_DIR \
						$CORE_PATH \
						$PROJECT_SAMPLE \
						$SM_TAG
				}

	##########################################################################
	### grabbing per sample indel only vcf files for on bait and on target ###
	##########################################################################

		# for each sample, make a vcf containing only passing mixed

			PASSING_MIXED_ON_BAIT_BY_SAMPLE ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A09_PASSING_MIXED_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A09_PASSING_MIXED_ON_BAIT_BY_SAMPLE_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A09_PASSING_MIXED_ON_BAIT_BY_SAMPLE.sh \
					$JAVA_1_8 \
					$GATK_DIR \
					$SAMPLE_REF_GENOME \
					$CORE_PATH \
					$PROJECT_SAMPLE \
					$SM_TAG
			}

		# for each sample, make a vcf containing only passing indels that fall within the on target bed file

			PASSING_MIXED_ON_TARGET_BY_SAMPLE ()
			{
				echo \
				qsub \
					-S /bin/bash \
					-cwd \
					-V \
					-q $QUEUE_LIST \
					-p $PRIORITY \
					-j y \
				-N K01A10_PASSING_MIXED_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
					-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/K01A10_PASSING_MIXED_ON_TARGET_BY_SAMPLE_$SAMPLE.log \
				-hold_jid K01_SELECT_VARIANTS_FOR_SAMPLE_$UNIQUE_ID_SM_TAG \
				$SCRIPT_DIR/K01A10_PASSING_MIXED_ON_TARGET_BY_SAMPLE.sh \
					$JAVA_1_8 \
					$GATK_DIR \
					$SAMPLE_REF_GENOME \
					$CORE_PATH \
					$PROJECT_SAMPLE \
					$SM_TAG \
					$TARGET_BED
			}

	# MAKE A QC REPORT FOR EACH SAMPLE
	# THIS IS CREATING A JOB_ID FOR A SAMPLE WHEN ALL OF THE BREAKOUTs PER SAMPLE IS DONE
	# THIS IS TO MITIGATE CREATING A HOLD ID THAT IS TOO LONG FOR GENERATING THE QC REPORT.
	# ALTHOUGH AT SOME POINT THIS STRING MIGHT END BEING TOO LONG AT SOME POINT.
	# SO QC REPORTS MIGHT HAVE TO END UP BEING DONE OUTSIDE OF THE PIPELINE FOR SOME BIG PROJECTS.
	# yucky, yuck...indenting creates a white space in the hold id which does not work so I have to do this hot mess.
	
		QC_REPORT_PREP ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
				-j y \
			-N Y"_"$BARCODE_2D \
				-o $CORE_PATH/$PROJECT_MS/LOGS/$SM_TAG/$SM_TAG"-QC_REPORT_PREP_QC.log" \
			-hold_jid K01A01_PASSING_VARIANTS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG,\
K01A02_PASSING_SNVS_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG,\
K01A03-1_CONCORDANCE_ON_TARGET_PER_SAMPLE_$UNIQUE_ID_SM_TAG,\
K01A04_PASSING_INDELS_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG,\
K01A05_PASSING_INDELS_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG,\
K01A06-1_TITV_ALL_$UNIQUE_ID_SM_TAG,\
K01A07-1_TITV_KNOWN_$UNIQUE_ID_SM_TAG,\
K01A08-1_TITV_NOVEL_$UNIQUE_ID_SM_TAG,\
K01A09_PASSING_MIXED_ON_BAIT_BY_SAMPLE_$UNIQUE_ID_SM_TAG,\
K01A10_PASSING_MIXED_ON_TARGET_BY_SAMPLE_$UNIQUE_ID_SM_TAG \
$SCRIPT_DIR/Y01_QC_REPORT_PREP.sh \
	$SAMTOOLS_DIR \
	$DATAMASH_DIR \
	$CORE_PATH \
	$PROJECT_SAMPLE \
	$SM_TAG \
	$PROJECT_MS \
	$PREFIX
		}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq )
do
	CREATE_SAMPLE_INFO_ARRAY
	MAKE_PROJ_DIR_TREE
	SELECT_PASSING_VARIANTS_PER_SAMPLE
	echo sleep 0.1s
	PASSING_VARIANTS_ON_TARGET_BY_SAMPLE
	echo sleep 0.1s
	PASSING_SNVS_ON_BAIT_BY_SAMPLE
	echo sleep 0.1s
	PASSING_SNVS_ON_TARGET_BY_SAMPLE
	echo sleep 0.1s
	CONCORDANCE_ON_TARGET_PER_SAMPLE
	echo sleep 0.1s
	PASSING_INDELS_ON_BAIT_BY_SAMPLE
	echo sleep 0.1s
	PASSING_INDELS_ON_TARGET_BY_SAMPLE
	echo sleep 0.1s
	PASSING_SNVS_TITV_ALL
	echo sleep 0.1s
	TITV_ALL
	echo sleep 0.1s
	PASSING_SNVS_TITV_KNOWN
	echo sleep 0.1s
	TITV_KNOWN
	echo sleep 0.1s
	PASSING_SNVS_TITV_NOVEL
	echo sleep 0.1s
	TITV_NOVEL
	echo sleep 0.1s
	PASSING_MIXED_ON_BAIT_BY_SAMPLE
	echo sleep 0.1s
	PASSING_MIXED_ON_TARGET_BY_SAMPLE
	echo sleep 0.1s
	QC_REPORT_PREP
	echo sleep 0.1s
done

##########################################################################
######################End of Functions####################################
##########################################################################

	# grab email addy

		SEND_TO=`cat $SCRIPT_DIR/../email_lists.txt`

	# Maybe I'll make this a function and throw it into a loop, but today is not that day.
	# I think that i will have to make this a look to handle multiple projects...maybe not
	# but again, today is not that day.

		awk 1 $SAMPLE_SHEET \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
			| awk 'BEGIN {FS=","; OFS="\t"} NR>1 {print $1,$8}' \
			| awk '{split($2,sm_tag,/[@-]/)} {print sm_tag[2]}' \
			| sort -k 1,1 \
			| uniq \
			| $DATAMASH_DIR/datamash \
				-s \
				collapse 1 \
			| awk 'gsub (/,/,",Y_",$1) \
				{print "qsub",\
					"-S /bin/bash",\
					"-cwd",\
					"-V",\
					"-q" , "'$QUEUE_LIST'",\
					"-p" , "'$PRIORITY'",\
					"-j y",\
					"-m","e",\
					"-M","'$SEND_TO'",\
				"-N" , "Y01-Y01-END_PROJECT_TASKS_" "'$PREFIX'",\
					"-o","'$CORE_PATH'" "/" "'$PROJECT_MS'" "/LOGS/Y01-Y01-" "'$PREFIX'" ".END_PROJECT_TASKS.log",\
				"-hold_jid" , "Y_" $1 ",A02-LAB_PREP_METRICS_" "'$PROJECT_MS'",\
				"'$SCRIPT_DIR'" "/Y01-Y01_END_PROJECT_TASKS.sh",\
					"'$CORE_PATH'",\
					"'$DATAMASH_DIR'",\
					"'$PROJECT_MS'",\
					"'$PREFIX'",\
					"'$SAMPLE_SHEET'",\
					"'$BEDTOOLS_DIR'",\
					"'$REF_SEQ_TRANSCRIPTS'" "\n" "sleep 0.1s"}'

# email when finished submitting

	SUBMITTER_ID=`whoami`

	PERSON_NAME=`getent passwd | awk 'BEGIN {FS=":"} $1=="'$SUBMITTER_ID'" {print $5}'`

	SCATTER_COUNT=`ls $CORE_PATH/$PROJECT_MS/TEMP/BED_FILE_SPLIT/BF*bed | wc -l`

	STUDY_COUNT=`awk '{print "basename",$1,".g.vcf.gz"}' $GVCF_LIST | bash | grep ^[0-9] | wc -l`

	HAPMAP_COUNT=`awk '{print "basename",$1,".g.vcf.gz"}' $GVCF_LIST | bash | grep -v ^[0-9] | wc -l`

	printf "$SAMPLE_SHEET\nhas finished submitting at\n`date`\nby `whoami`\nMULTI-SAMPLE VCF OUTPUT PROJECT IS:\n$PROJECT_MS\nVCF PREFIX IS:\n$PREFIX\nSCATTER IS $SCATTER_COUNT\n$TOTAL_SAMPLES samples called together\n$STUDY_COUNT study samples\n$HAPMAP_COUNT HapMap samples" \
		| mail -s "$PERSON_NAME has submitted HARD_FILTER_ALL_SUBMITTER.sh" \
			$SEND_TO
