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

	SAMTOOLS_DIR=$1 # needs to be a newer version to read a cram file. e.g. 1.6
	DATAMASH_DIR=$2

	CORE_PATH=$3
	PROJECT_SAMPLE=$4
	SM_TAG=$5
	PROJECT_MS=$6
	PREFIX=$7

# next script will cat everything together and add the header.

#############################################################
##### Grabbing the BAM header (for RG ID,PU,LB,etc) #########
#############################################################
#############################################################
##### THIS IS THE HEADER ####################################
##### "PROJECT_SAMPLE","SM_TAG","RG_PU","LIBRARY" ###########
##### "LIBRARY_PLATE","LIBRARY_WELL","LIBRARY_ROW" ##########
##### "LIBRARY_COLUMN","HYB_PLATE","HYB_WELL","HYB_ROW" #####
##### "HYB_COLUMN" ##########################################
#############################################################

	if [ -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS/RG_HEADER/$SM_TAG".RG_HEADER.txt" ]
		then
			cat $CORE_PATH/$PROJECT_SAMPLE/REPORTS/RG_HEADER/$SM_TAG".RG_HEADER.txt" \
					| $DATAMASH_DIR/datamash \
						-s \
						-g 1,2 \
						collapse 3 \
						unique 4 \
						unique 5 \
						unique 6 \
						unique 7 \
						unique 8 \
						unique 9 \
						unique 10 \
						unique 11 \
						unique 12 \
					| sed 's/,/;/g' \
					| $DATAMASH_DIR/datamash transpose \
			>| $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		elif [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS/RG_HEADER/$SM_TAG".RG_HEADER.txt" && -f $CORE_PATH/$PROJECT_SAMPLE/CRAM/$SM_TAG".cram" ]];
			then

				# grab field number for SM_TAG

					SM_FIELD=(`$SAMTOOLS_DIR/samtools view -H \
					$CORE_PATH/$PROJECT_SAMPLE/CRAM/$SM_TAG".cram" \
						| grep -m 1 ^@RG \
						| sed 's/\t/\n/g' \
						| cat -n \
						| sed 's/^ *//g' \
						| awk '$2~/^SM:/ {print $1}'`)

				# grab field number for PLATFORM_UNIT_TAG

					PU_FIELD=(`$SAMTOOLS_DIR/samtools view -H \
					$CORE_PATH/$PROJECT_SAMPLE/CRAM/$SM_TAG".cram" \
						| grep -m 1 ^@RG \
						| sed 's/\t/\n/g' \
						| cat -n \
						| sed 's/^ *//g' \
						| awk '$2~/^PU:/ {print $1}'`)

				# grab field number for LIBRARY_TAG

					LB_FIELD=(`$SAMTOOLS_DIR/samtools view -H \
					$CORE_PATH/$PROJECT_SAMPLE/CRAM/$SM_TAG".cram" \
						| grep -m 1 ^@RG \
						| sed 's/\t/\n/g' \
						| cat -n \
						| sed 's/^ *//g' \
						| awk '$2~/^LB:/ {print $1}'`)

				# Now grab the header and format
					# breaking out the library name into its parts is assuming that the format is...
					# fill in empty fields with NA thing (for loop in awk) is a lifesaver
					# https://unix.stackexchange.com/questions/53448/replacing-missing-value-blank-space-with-zero

					$SAMTOOLS_DIR/samtools view -H \
					$CORE_PATH/$PROJECT_SAMPLE/CRAM/$SM_TAG".cram" \
						| grep ^@RG \
						| awk \
							-v SM_FIELD="$SM_FIELD" \
							-v PU_FIELD="$PU_FIELD" \
							-v LB_FIELD="$LB_FIELD" \
							'BEGIN {OFS="\t"} {split($SM_FIELD,SMtag,":"); split($PU_FIELD,PU,":"); split($LB_FIELD,Library,":"); split(Library[2],Library_Unit,"_"); \
							print "'$PROJECT_SAMPLE'",SMtag[2],PU[2],Library[2],Library_Unit[1],Library_Unit[2],substr(Library_Unit[2],1,1),substr(Library_Unit[2],2,2),\
							Library_Unit[3],Library_Unit[4],substr(Library_Unit[4],1,1),substr(Library_Unit[4],2,2)}' \
						| awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1' \
						| $DATAMASH_DIR/datamash \
							-s \
							-g 1,2 \
							collapse 3 \
							unique 4 \
							unique 5 \
							unique 6 \
							unique 7 \
							unique 8 \
							unique 9 \
							unique 10 \
							unique 11 \
							unique 12 \
						| sed 's/,/;/g' \
						| $DATAMASH_DIR/datamash transpose \
					>| $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else
			echo -e "$PROJECT\t$SM_TAG\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" \
			| $DATAMASH_DIR/datamash transpose \
			>| $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

#################################################
##### GENDER CHECK FROM ANEUPLOIDY CHECK ########
#################################################
##### THIS IS THE HEADER ########################
##### X_AVG_DP,X_NORM_DP,Y_AVG_DP,Y_NORM_DP #####
#################################################

	awk 'BEGIN {OFS="\t"} $2=="X"&&$3=="whole" {print $6,$7} $2=="Y"&&$3=="whole" {print $6,$7}' \
	$CORE_PATH/$PROJECT_SAMPLE/REPORTS/ANEUPLOIDY_CHECK/$SM_TAG".chrom_count_report.txt" \
		| paste - - \
		| awk 'BEGIN {OFS="\t"} END {if ($1!~/[0-9]/) print "NaN","NaN","NaN","NaN"; else print $0}' \
		| $DATAMASH_DIR/datamash transpose \
	>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"

#################################
##### GRABBING CONCORDANCE. #####
#################################
##########################################################################
##### THIS IS THE HEADER #################################################
##### "COUNT_DISC_HOM","COUNT_CONC_HOM","PERCENT_CONC_HOM", ##############
##### "COUNT_DISC_HET","COUNT_CONC_HET","PERCENT_CONC_HET", ##############
##### "PERCENT_TOTAL_CONC","COUNT_HET_BEADCHIP","SENSITIVITY_2_HET" ######
##### "SNP_ARRAY" ########################################################
##########################################################################

	if [ -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS/CONCORDANCE_MS/$SM_TAG"_concordance.csv" ];
		then
			awk 1 $CORE_PATH/$PROJECT_SAMPLE/REPORTS/CONCORDANCE_MS/$SM_TAG"_concordance.csv" \
			| awk 'BEGIN {FS=",";OFS="\t"} NR>1 \
			{print $5,$6,$7,$2,$3,$4,$8,$9,$10,$11}' \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else
			echo -e "NaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN" \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

#####################################################################################################################################
##### VERIFY BAM ID #################################################################################################################
#####################################################################################################################################
##### THIS IS THE HEADER ############################################################################################################
##### "VERIFYBAM_FREEMIX","VERIFYBAM_#SNPS","VERIFYBAM_FREELK1","VERIFYBAM_FREELK0","VERIFYBAM_DIFF_LK0_LK1","VERIFYBAM_AVG_DP" #####
#####################################################################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//VERIFYBAMID/$SM_TAG".selfSM" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else
			awk 'BEGIN {OFS="\t"} NR>1 {print $7*100,$4,$8,$9,($9-$8),$6}' \
			$CORE_PATH/$PROJECT_SAMPLE/REPORTS//VERIFYBAMID/$SM_TAG".selfSM" \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

######################################################################################################
##### INSERT SIZE ####################################################################################
######################################################################################################
##### THIS IS THE HEADER #############################################################################
##### "MEDIAN_INSERT_SIZE","MEAN_INSERT_SIZE","STANDARD_DEVIATION_INSERT_SIZE","MAD_INSERT_SIZE" #####
######################################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS/INSERT_SIZE/METRICS/$SM_TAG".insert_size_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else
			awk 'BEGIN {OFS="\t"} NR==8 {print $1,$6,$7,$3}' \
			$CORE_PATH/$PROJECT_SAMPLE/REPORTS/INSERT_SIZE/METRICS/$SM_TAG".insert_size_metrics.txt" \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

#######################################################################################################
##### ALIGNMENT SUMMARY METRICS FOR READ 1 ############################################################
#######################################################################################################
##### THIS THE HEADER #################################################################################
##### "PCT_PF_READS_ALIGNED_R1","PF_HQ_ALIGNED_READS_R1","PF_HQ_ALIGNED_Q20_BASES_R1" #################
##### "PF_MISMATCH_RATE_R1","PF_HQ_ERROR_RATE_R1","PF_INDEL_RATE_R1" ##################################
##### "PCT_READS_ALIGNED_IN_PAIRS_R1","PCT_ADAPTER_R1" ################################################
#######################################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else
			awk 'BEGIN {OFS="\t"} NR==8 {if ($1=="UNPAIRED") print "0","0","0","0","0","0","0","0"; \
				else print $7*100,$9,$11,$13,$14,$15,$18*100,$24*100}' \
			$CORE_PATH/$PROJECT_SAMPLE/REPORTS//ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

#######################################################################################################
##### ALIGNMENT SUMMARY METRICS FOR READ 2 ############################################################
#######################################################################################################
##### THIS THE HEADER #################################################################################
##### "PCT_PF_READS_ALIGNED_R2","PF_HQ_ALIGNED_READS_R2","PF_HQ_ALIGNED_Q20_BASES_R2" #################
##### "PF_MISMATCH_RATE_R2","PF_HQ_ERROR_RATE_R2","PF_INDEL_RATE_R2" ##################################
##### "PCT_READS_ALIGNED_IN_PAIRS_R2","PCT_ADAPTER_R2" ################################################
#######################################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			awk 'BEGIN {OFS="\t"} NR==9 {if ($1=="") print "0","0","0","0","0","0","0","0" ; \
				else print $7*100,$9,$11,$13,$14,$15,$18*100,$24*100}' \
			$CORE_PATH/$PROJECT_SAMPLE/REPORTS//ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

################################################################################################
##### ALIGNMENT SUMMARY METRICS FOR PAIR #######################################################
################################################################################################
##### THIS THE HEADER ##########################################################################
##### "TOTAL_READS","RAW_GIGS","PCT_PF_READS_ALIGNED_PAIR" #####################################
##### "PF_MISMATCH_RATE_PAIR","PF_HQ_ERROR_RATE_PAIR","PF_INDEL_RATE_PAIR" #####################
##### "PCT_READS_ALIGNED_IN_PAIRS_PAIR","STRAND_BALANCE_PAIR","PCT_CHIMERAS_PAIR" ##############
##### "PF_HQ_ALIGNED_Q20_BASES_PAIR","MEAN_READ_LENGTH","PCT_PF_READS_IMPROPER_PAIRS_PAIR" #####
################################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			awk 'BEGIN {OFS="\t"} \
				NR==10 \
				{if ($1=="") print "0","0","0","0","0","0","0","0","0","0","0","0" ; \
				else print $2,($2*$16/1000000000),$7*100,$13,$14,$15,$18*100,$22,$23*100,$11,$16,$20*100}' \
			$CORE_PATH/$PROJECT_SAMPLE/REPORTS//ALIGNMENT_SUMMARY/$SM_TAG".alignment_summary_metrics.txt" \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

##################################################################################################################
##### MARK DUPLICATES REPORT #####################################################################################
##################################################################################################################
##### THIS IS THE HEADER #########################################################################################
##### "UNMAPPED_READS","READ_PAIR_OPTICAL_DUPLICATES","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE" #############
##### "SECONDARY_OR_SUPPLEMENTARY_READS","READ_PAIR_DUPLICATES","READ_PAIRS_EXAMINED","PAIRED_DUP_RATE" ##########
##### "UNPAIRED_READ_DUPLICATES","UNPAIRED_READS_EXAMINED","UNPAIRED_DUP_RATE","PERCENT_DUPLICATION_OPTICAL" #####
##################################################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS/PICARD_DUPLICATES/$SM_TAG"_MARK_DUPLICATES.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else
			MAX_RECORD=(`grep -n "^$" $CORE_PATH/$PROJECT_SAMPLE/REPORTS/PICARD_DUPLICATES/$SM_TAG"_MARK_DUPLICATES.txt" | awk 'BEGIN {FS=":"} NR==2 {print $1}'`)

			awk 'BEGIN {OFS="\t"} \
				NR>7&&NR<'$MAX_RECORD' \
				{if ($10!~/[0-9]/) print $5,$8,"NaN","NaN",$4,$7,$3,"NaN",$6,$2,"NaN" ; \
				else if ($10~/[0-9]/&&$2=="0") print $5,$8,$9*100,$10,$4,$7,$3,($7/$3),$6,$2,"NaN" ; \
				else print $5,$8,$9*100,$10,$4,$7,$3,($7/$3),$6,$2,($6/$2)}' \
			$CORE_PATH/$PROJECT_SAMPLE/REPORTS/PICARD_DUPLICATES/$SM_TAG"_MARK_DUPLICATES.txt" \
			| $DATAMASH_DIR/datamash sum 1 sum 2 mean 4 sum 5 sum 6 sum 7 sum 9 sum 10 \
			| awk 'BEGIN {OFS="\t"} \
				{if ($3!~/[0-9]/) print $1,$2,"NaN","NaN",$4,$5,$6,"NaN",$7,$8,"NaN","NaN" ; \
				else if ($3~/[0-9]/&&$1=="0") print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,"NaN",($2/$6)*100 ; \
				else if ($3~/[0-9]/&&$1!="0"&&$8=="0") print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,"NaN",($2/$6)*100 ; \
				else print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,($7/$8),($2/$6)*100}' \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

##########################################################################################################################################################################
##### HYBRIDIZATION SELECTION REPORT #####################################################################################################################################
##########################################################################################################################################################################
##### THIS IS THE HEADER #################################################################################################################################################
##### "GENOME_SIZE","BAIT_TERRITORY","TARGET_TERRITORY","PCT_PF_UQ_READS_ALIGNED" ########################################################################################
##### "PF_UQ_GIGS_ALIGNED","PCT_SELECTED_BASES","ON_BAIT_VS_SELECTED","MEAN_BAIT_COVERAGE","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE","MAX_TARGET_COVERAGE" #########
##### "ZERO_CVG_TARGETS_PCT","PCT_EXC_MAPQ","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET","FOLD_80_BASE_PENALTY" #######################################################################
##### "PCT_TARGET_BASES_1X","PCT_TARGET_BASES_2X","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X" #####
##### "PCT_TARGET_BASES_100X","HS_LIBRARY_SIZE","AT_DROPOUT","GC_DROPOUT","THEORETICAL_HET_SENSITIVITY","HET_SNP_Q","BAIT_SET","PCT_USABLE_BASES_ON_BAIT"} ###############
##########################################################################################################################################################################

	# this will take when there are no reads in the file...but i don't think that it will handle when there are reads, but none fall on target
	# the next time i that happens i'll fix this to handle it.

		if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt" ]]
			then
				echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
				| $DATAMASH_DIR/datamash transpose \
				>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
			else
				awk 'BEGIN {FS="\t";OFS="\t"} \
					NR==8 \
					{if ($12=="?"&&$44=="") print $2,$3,$4,"NaN",($14/1000000000),"NaN","NaN",$23,$24,$25,$29,"NaN","NaN","NaN","NaN","NaN",\
					$36,$37,$38,$39,$40,$41,$42,$43,"NaN",$51,$52,$53,$54,$1,"NaN" ; \
					else if ($12!="?"&&$44=="") print $2,$3,$4,$12*100,($14/1000000000),$19*100,$21,$23,$24,$25,$29*100,$31*100,\
					$32*100,$33*100,$34*100,$35,$36*100,$37*100,$38*100,$39*100,$40*100,$41*100,$42*100,$43*100,"NaN",$51,$52,$53,$54,$1,$26*100 ; \
					else print $2,$3,$4,$12*100,($14/1000000000),$19*100,$21,$23,$24,$25,$29*100,$31*100,$32*100,$33*100,$34*100,$35,\
					$36*100,$37*100,$38*100,$39*100,$40*100,$41*100,$42*100,$43*100,$44,$51,$52,$53,$54,$1,$26*100}' \
				$CORE_PATH/$PROJECT_SAMPLE/REPORTS//HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt" \
				| $DATAMASH_DIR/datamash transpose \
				>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		fi

	# this was supposed to be
	## if there are no reads, then print x
	## if there are no reads in the target area, the print y
	## else the data is fine and do as you intended.
	## however i no longer have anything to test this on...

		# awk 'BEGIN {FS="\t";OFS="\t"} \
		# NR==8 \
		# {if ($12=="?"&&$44=="") print $2,$3,$4,"NaN",($14/1000000000),"NaN","NaN",$22,$23,$24,$25,$29,"NaN","NaN","NaN","NaN",\
		# $36,$37,$38,$39,$40,$41,$42,$43,"NaN",$51,$52,$53,$54,$1,"NaN" ; \
		# else if ($12!="?") print $2,$3,$4,$12,($14/1000000000),$19,$21,$22,$23,$24,$25,$29,$31,$32,$33,$34,\
		# $36,$37,$38,$39,$40,$41,$42,$43,"NaN",$51,$52,$53,$54,$1,$26 ; \
		# else print $2,$3,$4,$12,($14/1000000000),$19,$21,$22,$23,$24,$25,$29,$31,$32,$33,$34,\
		# $36,$37,$38,$39,$40,$41,$42,$43,$44,$51,$52,$53,$54,$1,$26}' \
		# $CORE_PATH/$PROJECT_SAMPLE/REPORTS//HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt" \
		# | $DATAMASH_DIR/datamash transpose \
		# >> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"

##############################################
##### BAIT BIAS REPORT FOR Cref and Gref #####
##############################################
##### THIS IS THE HEADER #####################
##### Cref_Q,Gref_Q ##########################
##############################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//BAIT_BIAS/SUMMARY/$SM_TAG".bait_bias_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else
			grep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/REPORTS//BAIT_BIAS/SUMMARY/$SM_TAG".bait_bias_summary_metrics.txt" \
				| sed '/^$/d' \
				| awk 'BEGIN {OFS="\t"} $12=="Cref"||$12=="Gref" {print $5}' \
				| paste - - \
				| $DATAMASH_DIR/datamash collapse 1 collapse 2 \
				| sed 's/,/;/g' \
				| awk 'BEGIN {OFS="\t"} {print $0}' \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

############################################################
##### PRE-ADAPTER BIAS REPORT FOR Deamination and OxoG #####
############################################################
##### THIS IS THE HEADER ###################################
##### DEAMINATION_Q,OxoG_Q #################################
############################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//PRE_ADAPTER/SUMMARY/$SM_TAG".pre_adapter_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			grep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/REPORTS//PRE_ADAPTER/SUMMARY/$SM_TAG".pre_adapter_summary_metrics.txt" \
				| sed '/^$/d' \
				| awk 'BEGIN {OFS="\t"} $12=="Deamination"||$12=="OxoG" {print $5}' \
				| paste - - \
				| $DATAMASH_DIR/datamash collapse 1 collapse 2 \
				| sed 's/,/;/g' \
				| awk 'BEGIN {OFS="\t"} {print $0}' \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

###########################################################
##### BASE DISTRIBUTION REPORT AVERAGE FROM PER CYCLE #####
###########################################################
##### THIS IS THE HEADER ##################################
##### PCT_A,PCT_C,PCT_G,PCT_T,PCT_N #######################
###########################################################

	BASE_DISTIBUTION_BY_CYCLE_ROW_COUNT=(`wc -l $CORE_PATH/$PROJECT_SAMPLE/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/METRICS/$SM_TAG".base_distribution_by_cycle_metrics.txt"`)

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//BASE_DISTRIBUTION_BY_CYCLE/METRICS/$SM_TAG".base_distribution_by_cycle_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"

		elif [[ -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/METRICS/$SM_TAG".base_distribution_by_cycle_metrics.txt" && $BASE_DISTIBUTION_BY_CYCLE_ROW_COUNT -lt 8 ]]
			then
				echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
				| datamash transpose \
				>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"

		else
			sed '/^$/d' $CORE_PATH/$PROJECT_SAMPLE/REPORTS//BASE_DISTRIBUTION_BY_CYCLE/METRICS/$SM_TAG".base_distribution_by_cycle_metrics.txt" \
				| awk 'NR>6' \
				| $DATAMASH_DIR/datamash \
					mean 3 \
					mean 4 \
					mean 5 \
					mean 6 \
					mean 7 \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

############################################
##### BASE SUBSTITUTION RATE ###############
############################################
##### THIS IS THE HEADER ###################
##### PCT_A_to_C,PCT_A_to_G,PCT_A_to_T #####
##### PCT_C_to_A,PCT_C_to_G,PCT_C_to_T #####
############################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/REPORTS//ERROR_SUMMARY/$SM_TAG".error_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			sed '/^$/d' $CORE_PATH/$PROJECT_SAMPLE/REPORTS//ERROR_SUMMARY/$SM_TAG".error_summary_metrics.txt" \
				| awk 'NR>6 {print $6*100}' \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

############################################################
##### GENERATE COUNT PCT,IN DBSNP FOR ON BAIT SNVS #########
############################################################
##### THIS IS THE HEADER ###################################
##### "COUNT_SNV_ON_BAIT","PERCENT_SNV_ON_BAIT_SNP138" #####
############################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait_SNV.vcf.gz" ]]
		then
			echo -e NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			zgrep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait_SNV.vcf.gz" \
				| awk '{SNV_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
					END {if (SNV_COUNT!="") {print SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100} \
					else {print "0","NaN"}}' \
				| sed 's/ /\t/g' \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

#######################################################################################
##### GENERATE COUNT PCT,IN DBSNP FOR ON TARGET SNVS ##################################
#######################################################################################
##### THIS IS THE HEADER ##############################################################
##### "COUNT_SNV_ON_TARGET""\t""PERCENT_SNV_ON_TARGET_SNP138""\t""HET:HOM_TARGET" ##### 
#######################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_SNV.vcf.gz" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			zgrep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/SNV/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_SNV.vcf.gz" \
				| awk '{SNV_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} {HET_COUNT+=($10 ~ /^0\/1/)} {VAR_HOM+=($10 ~ /^1\/1/)} \
					END {if (SNV_COUNT!=""&&VAR_HOM!="0") print SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100,(HET_COUNT)/(VAR_HOM); \
					else if (SNV_COUNT!=""&&VAR_HOM=="0") print SNV_COUNT,(DBSNP_COUNT/SNV_COUNT)*100,"NaN"; \
					else print "0","NaN","NaN"}' \
				| sed 's/ /\t/g' \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

####################################################
##### GRABBING TI/TV ON UCSC CODING EXONS, ALL #####
####################################################
##### THIS IS THE HEADER ###########################
##### "ALL_TI_TV_COUNT","ALL_TI_TV_RATIO" ##########
####################################################

	awk 'BEGIN {OFS="\t"} END {if ($2!="") {print $2,$6} \
	else {print "0","NaN"}}' \
	$CORE_PATH/$PROJECT_SAMPLE/REPORTS/TI_TV_MS/$SM_TAG"_All_.titv.txt" \
	| $DATAMASH_DIR/datamash transpose \
	>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"

######################################################
##### GRABBING TI/TV ON UCSC CODING EXONS, KNOWN #####
######################################################
##### THIS IS THE HEADER #############################
##### "KNOWN_TI_TV_COUNT","KNOWN_TI_TV_RATIO" ########
######################################################

	awk 'BEGIN {OFS="\t"} END {if ($2!="") {print $2,$6} \
	else {print "0","NaN"}}' \
	$CORE_PATH/$PROJECT_SAMPLE/REPORTS/TI_TV_MS/$SM_TAG"_Known_.titv.txt" \
	| $DATAMASH_DIR/datamash transpose \
	>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"

######################################################
##### GRABBING TI/TV ON UCSC CODING EXONS, NOVEL #####
######################################################
##### THIS IS THE HEADER #############################
##### "NOVEL_TI_TV_COUNT","NOVEL_TI_TV_RATIO" ########
######################################################

	awk 'BEGIN {OFS="\t"} END {if ($2!="") {print $2,$6} \
	else {print "0","NaN"}}' \
	$CORE_PATH/$PROJECT_SAMPLE/REPORTS/TI_TV_MS/$SM_TAG"_Novel_.titv.txt" \
	| $DATAMASH_DIR/datamash transpose \
	>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"

#############################################################################################################################
##### INDEL METRICS ON BAIT #################################################################################################
#############################################################################################################################
##### THIS IS THE HEADER ####################################################################################################
##### "COUNT_ALL_INDEL_BAIT","ALL_INDEL_BAIT_PCT_SNP138","COUNT_BIALLELIC_INDEL_BAIT","BIALLELIC_INDEL_BAIT_PCT_SNP138" #####
#############################################################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/INDEL/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait_INDEL.vcf.gz" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			zgrep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/INDEL/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait_INDEL.vcf.gz" \
				| awk '{INDEL_COUNT++NR} \
				{INDEL_BIALLELIC+=($5!~",")} \
				{DBSNP_COUNT+=($3~"rs")} \
				{DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
				END {if (INDEL_BIALLELIC==""&&INDEL_COUNT=="") print "0","NaN","0","NaN"; \
				else if (INDEL_BIALLELIC==0&&INDEL_COUNT>=1) print INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,"0","NaN"; \
				else print INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100}' \
				| sed 's/ /\t/g' \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

#####################################################################################################################################
##### INDEL METRICS ON TARGET #######################################################################################################
#####################################################################################################################################
##### THIS IS THE HEADER ############################################################################################################
##### "COUNT_ALL_INDEL_TARGET","ALL_INDEL_TARGET_PCT_SNP138","COUNT_BIALLELIC_INDEL_TARGET","BIALLELIC_INDEL_TARGET_PCT_SNP138" #####
##### "BIALLELIC_ID_RATIO" ##########################################################################################################
#####################################################################################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/INDEL/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_INDEL.vcf.gz" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			zgrep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/INDEL/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_INDEL.vcf.gz" \
				| awk '{INDEL_COUNT++NR} \
					{INDEL_BIALLELIC+=($5!~",")} \
					{DBSNP_COUNT+=($3~"rs")} \
					{DBSNP_COUNT_BIALLELIC+=($3~"rs"&&$5!~",")} \
					{BIALLELIC_INSERTION+=(length($5)-length($4))>0&&$5!~","} \
					{BIALLELIC_DELETION+=(length($5)-length($4))<0&&$5!~","} \
					END {if (INDEL_COUNT=="") print "0","NaN","0","NaN","NaN"; \
					else if (INDEL_COUNT!=""&&INDEL_BIALLELIC==0) print INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,"0","NaN","NaN"; \
					else if (INDEL_COUNT!=""&&INDEL_BIALLELIC>=1&&BIALLELIC_DELETION==0) print INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100,"NaN"; \
					else print INDEL_COUNT,(DBSNP_COUNT/INDEL_COUNT)*100,INDEL_BIALLELIC,(DBSNP_COUNT_BIALLELIC/INDEL_BIALLELIC)*100,(BIALLELIC_INSERTION/BIALLELIC_DELETION)}' \
				| sed 's/ /\t/g' \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

#################################################################
##### BASIC METRICS FOR MIXED VARIANT TYPES ON BAIT #############
#################################################################
##### GENERATE COUNT PCT,IN DBSNP FOR ON BAIT MIXED VARIANT #####
##### THIS IS THE HEADER ########################################
##### "COUNT_MIXED_ON_BAIT","PERCENT_MIXED_ON_BAIT_SNP138"} #####
#################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/MIXED/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait_MIXED.vcf.gz" ]]
		then
			echo -e NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			zgrep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/MIXED/RELEASE/FILTERED_ON_BAIT/$SM_TAG"_MS_OnBait_MIXED.vcf.gz" \
				| awk '{MIXED_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
					END {if (MIXED_COUNT!="") print MIXED_COUNT,(DBSNP_COUNT/MIXED_COUNT)*100 ; \
					else print "0","NaN"}' \
				| sed 's/ /\t/g' \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

####################################################################
##### GENERATE COUNT PCT,IN DBSNP FOR ON TARGET MIXED VARIANT ######
####################################################################
##### THIS IS THE HEADER ###########################################
##### "COUNT_MIXED_ON_TARGET","PERCENT_MIXED_ON_TARGET_SNP138" #####
####################################################################

	if [[ ! -f $CORE_PATH/$PROJECT_SAMPLE/MIXED/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_MIXED.vcf.gz" ]]
		then
			echo -e NaN'\t'NaN \
			| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
		else

			zgrep -v "^#" $CORE_PATH/$PROJECT_SAMPLE/MIXED/RELEASE/FILTERED_ON_TARGET/$SM_TAG"_MS_OnTarget_MIXED.vcf.gz" \
				| awk '{MIXED_COUNT++NR} {DBSNP_COUNT+=($3~"rs")} \
					END {if (MIXED_COUNT!="") print MIXED_COUNT,(DBSNP_COUNT/MIXED_COUNT)*100 ; \
					else print "0","NaN"}' \
				| sed 's/ /\t/g' \
				| $DATAMASH_DIR/datamash transpose \
			>> $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt"
	fi

##########################################################################################################
##### tranpose from rows to list so these can be concatenated together for a project/batch QC report #####
##########################################################################################################

	cat $CORE_PATH/$PROJECT_MS/TEMP/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_TEMP.txt" \
	| $DATAMASH_DIR/datamash transpose \
	>| $CORE_PATH/$PROJECT_MS/REPORTS/QC_REPORT_PREP_$PREFIX/$SM_TAG".QC_REPORT_PREP.txt"
