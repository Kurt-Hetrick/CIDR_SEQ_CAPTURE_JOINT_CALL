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

	CORE_PATH=$1
	DATAMASH=$2

	PROJECT_MS=$3
	PREFIX=$4

	SAMPLE_SHEET=$5

##############################
##############################
##### MAKE THE QC REPORT #####
##############################
##############################

	TIMESTAMP=`date '+%F.%H-%M-%S'`

	SAMPLE_SHEET_NAME=`basename $SAMPLE_SHEET .csv`

	#######################################################################################
	##### combine all the individual qc reports for the PROJECT_MS and add the header #####
	#######################################################################################

		cat $CORE_PATH/$PROJECT_MS/REPORTS/QC_REPORT_PREP_$PREFIX/*.QC_REPORT_PREP.txt \
				| awk 'BEGIN {print "PROJECT",\
					"SM_TAG",\
					"RG_PU",\
					"LIBRARY",\
					"LIBRARY_PLATE",\
					"LIBRARY_WELL",\
					"LIBRARY_ROW",\
					"LIBRARY_COLUMN",\
					"HYB_PLATE",\
					"HYB_WELL",\
					"HYB_ROW",\
					"HYB_COLUMN",\
					"X_AVG_DP",\
					"X_NORM_DP",\
					"Y_AVG_DP",\
					"Y_NORM_DP",\
					"COUNT_DISC_HOM",\
					"COUNT_CONC_HOM",\
					"PERCENT_CONC_HOM",\
					"COUNT_DISC_HET",\
					"COUNT_CONC_HET",\
					"PERCENT_CONC_HET",\
					"PERCENT_TOTAL_CONC",\
					"COUNT_HET_BEADCHIP",\
					"SENSITIVITY_2_HET",\
					"SNP_ARRAY",\
					"VERIFYBAM_FREEMIX",\
					"VERIFYBAM_#SNPS",\
					"VERIFYBAM_FREELK1",\
					"VERIFYBAM_FREELK0",\
					"VERIFYBAM_DIFF_LK0_LK1",\
					"VERIFYBAM_AVG_DP",\
					"MEDIAN_INSERT_SIZE",\
					"MEAN_INSERT_SIZE",\
					"STANDARD_DEVIATION_INSERT_SIZE",\
					"MAD_INSERT_SIZE",\
					"PCT_PF_READS_ALIGNED_R1",\
					"PF_HQ_ALIGNED_READS_R1",\
					"PF_HQ_ALIGNED_Q20_BASES_R1",\
					"PF_MISMATCH_RATE_R1",\
					"PF_HQ_ERROR_RATE_R1",\
					"PF_INDEL_RATE_R1",\
					"PCT_READS_ALIGNED_IN_PAIRS_R1",\
					"PCT_ADAPTER_R1",\
					"PCT_PF_READS_ALIGNED_R2",\
					"PF_HQ_ALIGNED_READS_R2",\
					"PF_HQ_ALIGNED_Q20_BASES_R2",\
					"PF_MISMATCH_RATE_R2",\
					"PF_HQ_ERROR_RATE_R2",\
					"PF_INDEL_RATE_R2",\
					"PCT_READS_ALIGNED_IN_PAIRS_R2",\
					"PCT_ADAPTER_R2",\
					"TOTAL_READS",\
					"RAW_GIGS",\
					"PCT_PF_READS_ALIGNED_PAIR",\
					"PF_MISMATCH_RATE_PAIR",\
					"PF_HQ_ERROR_RATE_PAIR",\
					"PF_INDEL_RATE_PAIR",\
					"PCT_READS_ALIGNED_IN_PAIRS_PAIR",\
					"STRAND_BALANCE_PAIR",\
					"PCT_CHIMERAS_PAIR",\
					"PF_HQ_ALIGNED_Q20_BASES_PAIR",\
					"MEAN_READ_LENGTH",\
					"PCT_PF_READS_IMPROPER_PAIRS_PAIR",\
					"UNMAPPED_READS",\
					"READ_PAIR_OPTICAL_DUPLICATES",\
					"PERCENT_DUPLICATION",\
					"ESTIMATED_LIBRARY_SIZE",\
					"SECONDARY_OR_SUPPLEMENTARY_READS",\
					"READ_PAIR_DUPLICATES",\
					"READ_PAIRS_EXAMINED",\
					"PAIRED_DUP_RATE",\
					"UNPAIRED_READ_DUPLICATES",\
					"UNPAIRED_READS_EXAMINED",\
					"UNPAIRED_DUP_RATE",\
					"GENOME_SIZE",\
					"BAIT_TERRITORY",\
					"TARGET_TERRITORY",\
					"PCT_PF_UQ_READS_ALIGNED",\
					"PF_UQ_GIGS_ALIGNED",\
					"PCT_SELECTED_BASES",\
					"ON_BAIT_VS_SELECTED",\
					"MEAN_BAIT_COVERAGE",\
					"MEAN_TARGET_COVERAGE",\
					"MEDIAN_TARGET_COVERAGE",\
					"MAX_TARGET_COVERAGE",\
					"ZERO_CVG_TARGETS_PCT",\
					"PCT_EXC_MAPQ",\
					"PCT_EXC_BASEQ",\
					"PCT_EXC_OVERLAP",\
					"PCT_EXC_OFF_TARGET",\
					"PCT_TARGET_BASES_1X",\
					"PCT_TARGET_BASES_2X",\
					"PCT_TARGET_BASES_10X",\
					"PCT_TARGET_BASES_20X",\
					"PCT_TARGET_BASES_30X",\
					"PCT_TARGET_BASES_40X",\
					"PCT_TARGET_BASES_50X",\
					"PCT_TARGET_BASES_100X",\
					"HS_LIBRARY_SIZE",\
					"AT_DROPOUT",\
					"GC_DROPOUT",\
					"THEORETICAL_HET_SENSITIVITY",\
					"HET_SNP_Q",\
					"BAIT_SET",\
					"PCT_USABLE_BASES_ON_BAIT",\
					"Cref_Q",\
					"Gref_Q",\
					"DEAMINATION_Q",\
					"OxoG_Q",\
					"COUNT_SNV_ON_BAIT",\
					"PERCENT_SNV_ON_BAIT_SNP138",\
					"COUNT_SNV_ON_TARGET",\
					"PERCENT_SNV_ON_TARGET_SNP138",\
					"HET:HOM_TARGET",\
					"ALL_TI_TV_COUNT",\
					"ALL_TI_TV_RATIO",\
					"KNOWN_TI_TV_COUNT",\
					"KNOWN_TI_TV_RATIO",\
					"NOVEL_TI_TV_COUNT",\
					"NOVEL_TI_TV_RATIO",\
					"COUNT_ALL_INDEL_BAIT",\
					"ALL_INDEL_BAIT_PCT_SNP138",\
					"COUNT_BIALLELIC_INDEL_BAIT",\
					"BIALLELIC_INDEL_BAIT_PCT_SNP138",\
					"COUNT_ALL_INDEL_TARGET",\
					"ALL_INDEL_TARGET_PCT_SNP138",\
					"COUNT_BIALLELIC_INDEL_TARGET",\
					"BIALLELIC_INDEL_TARGET_PCT_SNP138",\
					"BIALLELIC_ID_RATIO",\
					"COUNT_MIXED_ON_BAIT",\
					"PERCENT_MIXED_ON_BAIT_SNP138",\
					"COUNT_MIXED_ON_TARGET",\
					"PERCENT_MIXED_ON_TARGET_SNP138"} \
					{print $0}' \
				| sed 's/ /,/g' \
				| sed 's/\t/,/g' \
		>| $CORE_PATH/$PROJECT_MS/TEMP/$PREFIX".QC_REPORT."$TIMESTAMP".TEMP.csv"

	#########################################################################
	##### Join with LAB QC PREP METRICS AND METADATA at the batch level #####
	#########################################################################

		join -t , -1 2 -2 1 \
		$CORE_PATH/$PROJECT_MS/TEMP/$PREFIX".QC_REPORT."$TIMESTAMP".TEMP.csv" \
		$CORE_PATH/$PROJECT_MS/REPORTS/LAB_PREP_REPORTS_MS/$SAMPLE_SHEET_NAME".LAB_PREP_METRICS.csv" \
		>| $CORE_PATH/$PROJECT_MS/REPORTS/QC_REPORTS/$SAMPLE_SHEET_NAME".QC_REPORT.csv"

######################################################################################################
######################################################################################################
##### MAKE ANEUPLOIDY AND PER CHROMOSOME VERIFYBAMID REPORTS FOR ALL SAMPLES IN THE SAMPLE SHEET #####
######################################################################################################
######################################################################################################

	# CREATE AND ARRAY FOR UNIQ PROJECT SAMPLE COMBINATIONS.

		CREATE_SAMPLE_ARRAY ()
		{
		SAMPLE_ARRAY=(`awk 1 $SAMPLE_SHEET \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} $8=="'$SM_TAG'" {print $1,$8}' \
			| sort \
			| uniq`)

			# $1
			PROJECT_SAMPLE=${SAMPLE_ARRAY[0]}

			# $8
			SM_TAG=${SAMPLE_ARRAY[1]}

		}

	# COMBINE ALL OF THE SAMPLE LEVEL ANEUPLOIDY AND PER CHROMOSOME VERIFYBAMID REPORTS INTO ONE FILE

		for SM_TAG in $(awk 'BEGIN {FS=","} $1=="'$PROJECT_MS'" {print $8}' $SAMPLE_SHEET \
			| sort \
			| uniq );
		do
				CREATE_SAMPLE_ARRAY

				######################################################
				#### Concatenate all aneuploidy reports together #####
				######################################################

					cat $CORE_PATH/$PROJECT_SAMPLE/REPORTS/ANEUPLOIDY_CHECK/$SM_TAG".chrom_count_report.txt" \
					>> $CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".ANEUPLOIDY_REPORT."$TIMESTAMP".txt"

				#######################################################################
				##### Concatenate all per chromosome verifybamID reports together #####
				#######################################################################

					cat $CORE_PATH/$PROJECT_SAMPLE/REPORTS/VERIFYBAMID_CHR/$SM_TAG".VERIFYBAMID.PER_CHR.txt " \
					>> $CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".PER_CHR_VERIFYBAMID."$TIMESTAMP".txt"		

		done

# FORMAT THE ANEUPLOIDY CHECK REPORT

	( cat $CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".ANEUPLOIDY_REPORT."$TIMESTAMP".txt" | grep "^SM_TAG" | uniq ; \
		cat $CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".ANEUPLOIDY_REPORT."$TIMESTAMP".txt" | grep -v "SM_TAG" ) \
		| sed 's/\t/,/g' \
		>| $CORE_PATH/$PROJECT_MS/REPORTS/QC_REPORTS/$PREFIX".ANEUPLOIDY_CHECK."$TIMESTAMP".csv"

# FORMAT THE PER CHROMOSOME VERIFYBAMID REPORT

	( cat $CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".PER_CHR_VERIFYBAMID."$TIMESTAMP".txt" | grep "^#" | uniq ; \
			cat $CORE_PATH/$PROJECT_MS/TEMP/$SAMPLE_SHEET_NAME".PER_CHR_VERIFYBAMID."$TIMESTAMP".txt" | grep -v "^#" ) \
			| sed 's/\t/,/g' \
			>| $CORE_PATH/$PROJECT_MS/REPORTS/QC_REPORTS/$PREFIX".PER_CHR_VERIFYBAMID."$TIMESTAMP".csv"

###################################################
#### Clean up the Wall Clock minutes tracker. #####
###################################################

	awk 'BEGIN {FS=",";OFS=","} $1~/^[A-Z 0-9]/&&$2!=""&&$3!=""&&$4!=""&&$5!=""&&$6!=""&&$7==""&&$5!~/A-Z/&&$6!~/A-Z/ \
	{print $1,$2,$3,$4,$5,$6,($6-$5)/60,strftime("%F.%H-%M-%S",$5),strftime("%F.%H-%M-%S",$6)}' \
	$CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv" \
		| awk 'BEGIN {print "SAMPLE_GROUP,TASK_GROUP,TASK,HOST,EPOCH_START,EPOCH_END,WC_MIN,TIMESTAMP_START,TIMESTAMP_END"} \
		{print $0}' \
	>|$CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.FIXED.csv"

#############################################################
##### Summarize Wall Clock times ############################
##### This is probably garbage. I'll look at this later #####
#############################################################

# sed 's/,/\t/g' $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".WALL.CLOCK.TIMES.csv" \
# | sort -k 1,1 -k 2,2 -k 3,3 \
# | awk 'BEGIN {OFS="\t"} {print $0,($6-$5),($6-$5)/60,($6-$5)/3600}' \
# | $DATAMASH/datamash -s -g 1,2 max 7 max 8 max 9 | tee $CORE_PATH/$PROJECT_MS/TEMP/WALL.CLOCK.TIMES.BY.GROUP.txt \
# | $DATAMASH/datamash -g 1 sum 3 sum 4 sum 5 \
# | awk 'BEGIN {print "SAMPLE_PROJECT_MS","WALL_CLOCK_SECONDS","WALL_CLOCK_MINUTES","WALL_CLOCK_HOURS"} {print $0}' \
# | sed -r 's/[[:space:]]+/,/g' \
# >| $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".WALL.CLOCK.TIMES.BY_SAMPLE.csv"

# sed 's/\t/,/g' $CORE_PATH/$PROJECT_MS/TEMP/WALL.CLOCK.TIMES.BY.GROUP.txt \
# | awk 'BEGIN {print "SAMPLE_PROJECT_MS","TASK_GROUP","WALL_CLOCK_SECONDS","WALL_CLOCK_MINUTES","WALL_CLOCK_HOURS"} {print $0}' \
# | sed -r 's/[[:space:]]+/,/g' \
# >| $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".WALL.CLOCK.TIMES.BY_SAMPLE_GROUP.csv"

# echo PROJECT_MS finished at `date` >> $CORE_PATH/$PROJECT_MS/REPORTS/PROJECT_MS_START_END_TIMESTAMP.txt

# todo: oxidation report?
