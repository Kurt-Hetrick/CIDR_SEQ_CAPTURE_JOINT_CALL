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
	PARALLEL_DIR=$3

	CORE_PATH=$4
	PROJECT_MS=$5
	PREFIX=$6
	SAMPLE_SHEET=$7

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

echo
echo RETRIEVING READ GROUP HEADERS: `date`
echo

	GRAB_READ_GROUP_HEADER ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6
		SAMTOOLS_DIR=$7

		if [ -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/RG_HEADER/${SM_TAG}.RG_HEADER.txt ]
		then
			cat ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/RG_HEADER/${SM_TAG}.RG_HEADER.txt \
					| ${DATAMASH_DIR}/datamash \
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
					| ${DATAMASH_DIR}/datamash \
						transpose \
			>| ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}.QC_REPORT_TEMP.txt
		elif [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/RG_HEADER/${SM_TAG}.RG_HEADER.txt && -f ${CORE_PATH}/${PROJECT_SAMPLE}/CRAM/${SM_TAG}.cram ]];
			then

				# grab field number for SM_TAG

					SM_FIELD=(`${SAMTOOLS_DIR}/samtools \
						view -H \
					${CORE_PATH}/${PROJECT_SAMPLE}/CRAM/${SM_TAG}.cram \
						| grep -m 1 ^@RG \
						| sed 's/\t/\n/g' \
						| cat -n \
						| sed 's/^ *//g' \
						| awk '$2~/^SM:/ {print $1}'`)

				# grab field number for PLATFORM_UNIT_TAG

					PU_FIELD=(`${SAMTOOLS_DIR}/samtools \
						view -H \
					${CORE_PATH}/${PROJECT_SAMPLE}/CRAM/${SM_TAG}.cram \
						| grep -m 1 ^@RG \
						| sed 's/\t/\n/g' \
						| cat -n \
						| sed 's/^ *//g' \
						| awk '$2~/^PU:/ {print $1}'`)

				# grab field number for LIBRARY_TAG

					LB_FIELD=(`${SAMTOOLS_DIR}/samtools \
						view -H \
					${CORE_PATH}/${PROJECT_SAMPLE}/CRAM/${SM_TAG}.cram \
						| grep -m 1 ^@RG \
						| sed 's/\t/\n/g' \
						| cat -n \
						| sed 's/^ *//g' \
						| awk '$2~/^LB:/ {print $1}'`)

				# Now grab the header and format
					# breaking out the library name into its parts is assuming that the format is...
					# fill in empty fields with NA thing (for loop in awk) is a lifesaver
					# https://unix.stackexchange.com/questions/53448/replacing-missing-value-blank-space-with-zero

					${SAMTOOLS_DIR}/samtools \
						view -H \
					${CORE_PATH}/${PROJECT_SAMPLE}/CRAM/${SM_TAG}.cram \
						| grep ^@RG \
						| awk \
							-v SM_FIELD="$SM_FIELD" \
							-v PU_FIELD="$PU_FIELD" \
							-v LB_FIELD="$LB_FIELD" \
							'BEGIN {OFS="\t"} {split($SM_FIELD,SMtag,":"); split($PU_FIELD,PU,":"); split($LB_FIELD,Library,":"); split(Library[2],Library_Unit,"_"); \
							print "'${PROJECT_SAMPLE}'",SMtag[2],PU[2],Library[2],Library_Unit[1],Library_Unit[2],substr(Library_Unit[2],1,1),substr(Library_Unit[2],2,2),\
							Library_Unit[3],Library_Unit[4],substr(Library_Unit[4],1,1),substr(Library_Unit[4],2,2)}' \
						| awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1' \
						| ${DATAMASH_DIR}/datamash \
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
						| ${DATAMASH_DIR}/datamash \
							transpose \
					>| ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}.QC_REPORT_TEMP.txt
		else
			echo -e "${PROJECT_SAMPLE}\t${SM_TAG}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" \
			| ${DATAMASH_DIR}/datamash \
				transpose \
			>| ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}.QC_REPORT_TEMP.txt
		fi
	}

		export -f GRAB_READ_GROUP_HEADER

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_READ_GROUP_HEADER \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX \
			$SAMTOOLS_DIR

#################################################
##### GENDER CHECK FROM ANEUPLOIDY CHECK ########
#################################################
##### THIS IS THE HEADER ########################
##### X_AVG_DP,X_NORM_DP,Y_AVG_DP,Y_NORM_DP #####
#################################################

echo RETRIEVING GENDER CHECK METRICS FROM ANEUPLOIDY CHECK METRICS: `date`
echo

	GRAB_GENDER_CHECK_FROM_ANEUPLOIDY_CHECK ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		awk 'BEGIN {OFS="\t"} $2=="X"&&$3=="whole" {print $6,$7} $2=="Y"&&$3=="whole" {print $6,$7}' \
		${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/ANEUPLOIDY_CHECK/${SM_TAG}".chrom_count_report.txt" \
			| paste - - \
			| awk 'BEGIN {OFS="\t"} END {if ($1!~/[0-9]/) print "NaN","NaN","NaN","NaN"; else print $0}' \
			| ${DATAMASH_DIR}/datamash transpose \
		>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
	}

		export -f GRAB_GENDER_CHECK_FROM_ANEUPLOIDY_CHECK

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_GENDER_CHECK_FROM_ANEUPLOIDY_CHECK \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

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

echo RETRIEVING CONCORDANCE METRICS: `date`
echo

	GRAB_CONCORDANCE ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [ -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/CONCORDANCE_MS/${SM_TAG}"_concordance.csv" ];
		then
			awk 1 ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/CONCORDANCE_MS/${SM_TAG}"_concordance.csv" \
			| awk 'BEGIN {FS=",";OFS="\t"} NR>1 \
			{print $5,$6,$7,$2,$3,$4,$8,$9,$10,$11}' \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			echo -e "NaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN" \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_CONCORDANCE

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_CONCORDANCE \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

#########################################################################################
##### VERIFY BAM ID #####################################################################
#########################################################################################
##### THIS IS THE HEADER ################################################################
##### "VERIFYBAM_FREEMIX","VERIFYBAM_#SNPS","VERIFYBAM_FREELK1","VERIFYBAM_FREELK0" #####
##### "VERIFYBAM_DIFF_LK0_LK1","VERIFYBAM_AVG_DP" #######################################
#########################################################################################

echo RETRIEVING VERIFYBAMID METRICS: `date`
echo

	GRAB_VERIFYBAMID ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//VERIFYBAMID/${SM_TAG}".selfSM" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			awk 'BEGIN {OFS="\t"} NR>1 {print $7*100,$4,$8,$9,($9-$8),$6}' \
			${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//VERIFYBAMID/${SM_TAG}".selfSM" \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_VERIFYBAMID

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_VERIFYBAMID \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

######################################################################################################
##### INSERT SIZE ####################################################################################
######################################################################################################
##### THIS IS THE HEADER #############################################################################
##### "MEDIAN_INSERT_SIZE","MEAN_INSERT_SIZE","STANDARD_DEVIATION_INSERT_SIZE","MAD_INSERT_SIZE" #####
######################################################################################################

echo RETRIEVING INSERT SIZE METRICS: `date`
echo

	GRAB_INSERT_SIZE ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/INSERT_SIZE/METRICS/${SM_TAG}".insert_size_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			awk 'BEGIN {OFS="\t"} NR==8 {print $1,$6,$7,$3}' \
			${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/INSERT_SIZE/METRICS/${SM_TAG}".insert_size_metrics.txt" \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_INSERT_SIZE

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_INSERT_SIZE \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

#######################################################################################################
##### ALIGNMENT SUMMARY METRICS FOR READ 1 ############################################################
#######################################################################################################
##### THIS THE HEADER #################################################################################
##### "PCT_PF_READS_ALIGNED_R1","PF_HQ_ALIGNED_READS_R1","PF_HQ_ALIGNED_Q20_BASES_R1" #################
##### "PF_MISMATCH_RATE_R1","PF_HQ_ERROR_RATE_R1","PF_INDEL_RATE_R1" ##################################
##### "PCT_READS_ALIGNED_IN_PAIRS_R1","PCT_ADAPTER_R1" ################################################
#######################################################################################################

echo RETRIEVING ALIGNMENT SUMMARY METRICS FOR READ ONE: `date`
echo

	GRAB_ALIGNMENT_SUMMARY_METRICS_READ_ONE ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//ALIGNMENT_SUMMARY/${SM_TAG}".alignment_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			awk 'BEGIN {OFS="\t"} NR==8 {if ($1=="UNPAIRED") print "0","0","0","0","0","0","0","0"; \
				else print $7*100,$9,$11,$13,$14,$15,$18*100,$24*100}' \
			${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//ALIGNMENT_SUMMARY/${SM_TAG}".alignment_summary_metrics.txt" \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_ALIGNMENT_SUMMARY_METRICS_READ_ONE

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_ALIGNMENT_SUMMARY_METRICS_READ_ONE \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

#######################################################################################################
##### ALIGNMENT SUMMARY METRICS FOR READ 2 ############################################################
#######################################################################################################
##### THIS THE HEADER #################################################################################
##### "PCT_PF_READS_ALIGNED_R2","PF_HQ_ALIGNED_READS_R2","PF_HQ_ALIGNED_Q20_BASES_R2" #################
##### "PF_MISMATCH_RATE_R2","PF_HQ_ERROR_RATE_R2","PF_INDEL_RATE_R2" ##################################
##### "PCT_READS_ALIGNED_IN_PAIRS_R2","PCT_ADAPTER_R2" ################################################
#######################################################################################################

echo RETRIEVING READ ALIGNMENT SUMMARY METRICS FOR READ 2: `date`
echo

	GRAB_ALIGNMENT_SUMMARY_METRICS_READ_TWO ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//ALIGNMENT_SUMMARY/${SM_TAG}".alignment_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else

			awk 'BEGIN {OFS="\t"} NR==9 {if ($1=="") print "0","0","0","0","0","0","0","0" ; \
				else print $7*100,$9,$11,$13,$14,$15,$18*100,$24*100}' \
			${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//ALIGNMENT_SUMMARY/${SM_TAG}".alignment_summary_metrics.txt" \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_ALIGNMENT_SUMMARY_METRICS_READ_TWO

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_ALIGNMENT_SUMMARY_METRICS_READ_TWO \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

################################################################################################
##### ALIGNMENT SUMMARY METRICS FOR PAIR #######################################################
################################################################################################
##### THIS THE HEADER ##########################################################################
##### "TOTAL_READS","RAW_GIGS","PCT_PF_READS_ALIGNED_PAIR" #####################################
##### "PF_MISMATCH_RATE_PAIR","PF_HQ_ERROR_RATE_PAIR","PF_INDEL_RATE_PAIR" #####################
##### "PCT_READS_ALIGNED_IN_PAIRS_PAIR","STRAND_BALANCE_PAIR","PCT_CHIMERAS_PAIR" ##############
##### "PF_HQ_ALIGNED_Q20_BASES_PAIR","MEAN_READ_LENGTH","PCT_PF_READS_IMPROPER_PAIRS_PAIR" #####
################################################################################################

echo RETRIEVING READ ALIGNMENT SUMMARY METRICS FOR READ PAIRS: `date`
echo

	GRAB_ALIGNMENT_SUMMARY_METRICS_BY_PAIR ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//ALIGNMENT_SUMMARY/${SM_TAG}".alignment_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			awk 'BEGIN {OFS="\t"} \
				NR==10 \
				{if ($1=="") print "0","0","0","0","0","0","0","0","0","0","0","0" ; \
				else print $2,($2*$16/1000000000),$7*100,$13,$14,$15,$18*100,$22,$23*100,$11,$16,$20*100}' \
			${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//ALIGNMENT_SUMMARY/${SM_TAG}".alignment_summary_metrics.txt" \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_ALIGNMENT_SUMMARY_METRICS_BY_PAIR

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_ALIGNMENT_SUMMARY_METRICS_BY_PAIR \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

##################################
##### MARK DUPLICATES REPORT #####
##################################
##### THIS IS THE HEADER ####################################################################################
##### "UNMAPPED_READS","READ_PAIR_OPTICAL_DUPLICATES","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE" ########
##### "SECONDARY_OR_SUPPLEMENTARY_READS","READ_PAIR_DUPLICATES","READ_PAIRS_EXAMINED","PAIRED_DUP_RATE" #####
##### "UNPAIRED_READ_DUPLICATES","UNPAIRED_READS_EXAMINED","UNPAIRED_DUP_RATE" ##############################
##### "PERCENT_DUPLICATION_OPTICAL" #########################################################################
#############################################################################################################

echo RETRIEVING MARK DUPLICATES METRICS: `date`
echo

	GRAB_MARK_DUPLICATES ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/PICARD_DUPLICATES/${SM_TAG}"_MARK_DUPLICATES.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			MAX_RECORD=(`grep -n "^$" ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/PICARD_DUPLICATES/${SM_TAG}"_MARK_DUPLICATES.txt" | awk 'BEGIN {FS=":"} NR==2 {print $1}'`)

			awk 'BEGIN {OFS="\t"} \
				NR>7&&NR<'$MAX_RECORD' \
				{if ($10!~/[0-9]/) print $5,$8,"NaN","NaN",$4,$7,$3,"NaN",$6,$2,"NaN" ; \
				else if ($10~/[0-9]/&&$2=="0") print $5,$8,$9*100,$10,$4,$7,$3,($7/$3),$6,$2,"NaN" ; \
				else print $5,$8,$9*100,$10,$4,$7,$3,($7/$3),$6,$2,($6/$2)}' \
			${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/PICARD_DUPLICATES/${SM_TAG}"_MARK_DUPLICATES.txt" \
			| ${DATAMASH_DIR}/datamash sum 1 sum 2 mean 4 sum 5 sum 6 sum 7 sum 9 sum 10 \
			| awk 'BEGIN {OFS="\t"} \
				{if ($3!~/[0-9]/) print $1,$2,"NaN","NaN",$4,$5,$6,"NaN",$7,$8,"NaN","NaN" ; \
				else if ($3~/[0-9]/&&$1=="0") print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,"NaN",($2/$6)*100 ; \
				else if ($3~/[0-9]/&&$1!="0"&&$8=="0") print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,"NaN",($2/$6)*100 ; \
				else print $1,$2,(($7+($5*2))/($8+($6*2)))*100,$3,$4,$5,$6,($5/$6),$7,$8,($7/$8),($2/$6)*100}' \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_MARK_DUPLICATES

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_MARK_DUPLICATES \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

##########################################
##### HYBRIDIZATION SELECTION REPORT #####
##########################################
##### THIS IS THE HEADER ###################################################################################
##### "GENOME_SIZE","BAIT_TERRITORY","TARGET_TERRITORY","PCT_PF_UQ_READS_ALIGNED","PF_UQ_GIGS_ALIGNED" #####
##### "PCT_SELECTED_BASES","ON_BAIT_VS_SELECTED","MEAN_BAIT_COVERAGE","MEAN_TARGET_COVERAGE" ###############
##### "MEDIAN_TARGET_COVERAGE","MAX_TARGET_COVERAGE","ZERO_CVG_TARGETS_PCT","PCT_EXC_MAPQ" #################
##### "PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_OFF_TARGET","FOLD_80_BASE_PENALTY" ########################
##### "PCT_TARGET_BASES_1X","PCT_TARGET_BASES_2X","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X" ############
##### "PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X" #################################
##### "PCT_TARGET_BASES_100X","HS_LIBRARY_SIZE","AT_DROPOUT","GC_DROPOUT" ##################################
##### "THEORETICAL_HET_SENSITIVITY","HET_SNP_Q","BAIT_SET","PCT_USABLE_BASES_ON_BAIT" ######################
############################################################################################################

echo RETRIEVING HYBRIDIZATION METRICS: `date`
echo

	GRAB_HYB_SELECTION ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		# this will take when there are no reads in the file...but i don't think that it will handle when there are reads, but none fall on target
		# the next time i that happens i'll fix this to handle it.

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//HYB_SELECTION/${SM_TAG}"_hybridization_selection_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			awk 'BEGIN {FS="\t";OFS="\t"} \
				NR==8 \
				{if ($12=="?"&&$44=="") print $2,$3,$4,"NaN",($14/1000000000),"NaN","NaN",$23,$24,$25,$29,"NaN","NaN","NaN","NaN","NaN",\
				$36,$37,$38,$39,$40,$41,$42,$43,"NaN",$51,$52,$53,$54,$1,"NaN" ; \
				else if ($12!="?"&&$44=="") print $2,$3,$4,$12*100,($14/1000000000),$19*100,$21,$23,$24,$25,$29*100,$31*100,\
				$32*100,$33*100,$34*100,$35,$36*100,$37*100,$38*100,$39*100,$40*100,$41*100,$42*100,$43*100,"NaN",$51,$52,$53,$54,$1,$26*100 ; \
				else print $2,$3,$4,$12*100,($14/1000000000),$19*100,$21,$23,$24,$25,$29*100,$31*100,$32*100,$33*100,$34*100,$35,\
				$36*100,$37*100,$38*100,$39*100,$40*100,$41*100,$42*100,$43*100,$44,$51,$52,$53,$54,$1,$26*100}' \
			${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//HYB_SELECTION/${SM_TAG}"_hybridization_selection_metrics.txt" \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_HYB_SELECTION

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_HYB_SELECTION \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

##############################################
##### BAIT BIAS REPORT FOR Cref and Gref #####
##############################################
##### THIS IS THE HEADER #####################
##### Cref_Q,Gref_Q ##########################
##############################################

echo RETRIEVING BAIT BIAS METRICS: `date`
echo

	GRAB_BAIT_BIAS ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//BAIT_BIAS/SUMMARY/${SM_TAG}".bait_bias_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			grep -v "^#" ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//BAIT_BIAS/SUMMARY/${SM_TAG}".bait_bias_summary_metrics.txt" \
				| sed '/^$/d' \
				| awk 'BEGIN {OFS="\t"} $12=="Cref"||$12=="Gref" {print $5}' \
				| paste - - \
				| ${DATAMASH_DIR}/datamash collapse 1 collapse 2 \
				| sed 's/,/;/g' \
				| awk 'BEGIN {OFS="\t"} {print $0}' \
				| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_BAIT_BIAS

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_BAIT_BIAS \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

############################################################
##### PRE-ADAPTER BIAS REPORT FOR Deamination and OxoG #####
############################################################
##### THIS IS THE HEADER ###################################
##### DEAMINATION_Q,OxoG_Q #################################
############################################################

echo RETRIEVING PRE ADAPTER BIAS METRICS: `date`
echo

	GRAB_PRE_ADAPTER_BIAS ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//PRE_ADAPTER/SUMMARY/${SM_TAG}".pre_adapter_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else
			grep -v "^#" ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//PRE_ADAPTER/SUMMARY/${SM_TAG}".pre_adapter_summary_metrics.txt" \
				| sed '/^$/d' \
				| awk 'BEGIN {OFS="\t"} $12=="Deamination"||$12=="OxoG" {print $5}' \
				| paste - - \
				| ${DATAMASH_DIR}/datamash collapse 1 collapse 2 \
				| sed 's/,/;/g' \
				| awk 'BEGIN {OFS="\t"} {print $0}' \
				| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_PRE_ADAPTER_BIAS

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_PRE_ADAPTER_BIAS \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

###########################################################
##### BASE DISTRIBUTION REPORT AVERAGE FROM PER CYCLE #####
###########################################################
##### THIS IS THE HEADER ##################################
##### PCT_A,PCT_C,PCT_G,PCT_T,PCT_N #######################
###########################################################

echo RETRIEVING BASE DISTRIBUTION AVERAGE PER CYCLE METRICS: `date`
echo

	GRAB_BASE_DISTRIBUTION_AVERAGE_PER_CYCLE ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		BASE_DISTIBUTION_BY_CYCLE_ROW_COUNT=(`wc -l ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/METRICS/${SM_TAG}".base_distribution_by_cycle_metrics.txt"`)

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//BASE_DISTRIBUTION_BY_CYCLE/METRICS/${SM_TAG}".base_distribution_by_cycle_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"

		elif [[ -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/METRICS/${SM_TAG}".base_distribution_by_cycle_metrics.txt" && $BASE_DISTIBUTION_BY_CYCLE_ROW_COUNT -lt 8 ]]
			then
				echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
				| datamash transpose \
				>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"

		else
			sed '/^$/d' ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//BASE_DISTRIBUTION_BY_CYCLE/METRICS/${SM_TAG}".base_distribution_by_cycle_metrics.txt" \
				| awk 'NR>6' \
				| ${DATAMASH_DIR}/datamash \
					mean 3 \
					mean 4 \
					mean 5 \
					mean 6 \
					mean 7 \
				| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_BASE_DISTRIBUTION_AVERAGE_PER_CYCLE

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_BASE_DISTRIBUTION_AVERAGE_PER_CYCLE \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

############################################
##### BASE SUBSTITUTION RATE ###############
############################################
##### THIS IS THE HEADER ###################
##### PCT_A_to_C,PCT_A_to_G,PCT_A_to_T #####
##### PCT_C_to_A,PCT_C_to_G,PCT_C_to_T #####
############################################

echo RETRIEVING BASE SUBSTITUTION RATE METRICS: `date`
echo

	GRAB_BASE_SUBSTITUTION_RATE ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		if [[ ! -f ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//ERROR_SUMMARY/${SM_TAG}".error_summary_metrics.txt" ]]
		then
			echo -e NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN'\t'NaN \
			| ${DATAMASH_DIR}/datamash transpose \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		else

			sed '/^$/d' ${CORE_PATH}/${PROJECT_SAMPLE}/REPORTS//ERROR_SUMMARY/${SM_TAG}".error_summary_metrics.txt" \
				| awk 'NR>6 {print $6*100}' \
			>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt"
		fi
	}

		export -f GRAB_BASE_SUBSTITUTION_RATE

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_BASE_SUBSTITUTION_RATE \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

#############################################
##### VCF METRICS FOR BAIT BED FILE #########
#############################################
##### THIS IS THE HEADER ###########################################################################
##### "BAIT_HET_HOMVAR_RATIO","BAIT_PCT_GQ0_VARIANT","BAIT_TOTAL_GQ0_VARIANT" ######################
##### "BAIT_TOTAL_HET_DEPTH_SNV","BAIT_TOTAL_SNV","BAIT_NUM_IN_DBSNP_138_SNV","BAIT_NOVEL_SNV" #####
##### "BAIT_FILTERED_SNV","BAIT_PCT_DBSNP_138_SNV","BAIT_TOTAL_INDEL","BAIT_NOVEL_INDEL" ###########
##### "BAIT_FILTERED_INDEL","BAIT_PCT_DBSNP_138_INDEL","BAIT_NUM_IN_DBSNP_138_INDEL" ###############
##### "BAIT_DBSNP_138_INS_DEL_RATIO","BAIT_NOVEL_INS_DEL_RATIO","BAIT_TOTAL_MULTIALLELIC_SNV" ######
##### "BAIT_NUM_IN_DBSNP_138_MULTIALLELIC_SNV","BAIT_TOTAL_COMPLEX_INDEL" ##########################
##### "BAIT_NUM_IN_DBSNP_138_COMPLEX_INDEL","BAIT_SNP_REFERENCE_BIAS","BAIT_NUM_SINGLETONS" ########
####################################################################################################

echo RETRIEVING ON BAIT VCF METRICS: `date`
echo

	GRAB_VCF_METRICS_BAIT ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		awk 'BEGIN {OFS="\t"} \
				$1=="'${SM_TAG}'" \
				{print $2,$3*100,$4,$5,$6,$7,$8,$9,$10*100,$13,$14,$15,$16*100,\
					$17,$18,$19,$20,$21,$22,$23,$24,$25}' \
			${CORE_PATH}/${PROJECT_MS}/REPORTS/VCF_METRICS/MULTI_SAMPLE/${PREFIX}_BAIT.variant_calling_detail_metrics.txt \
			| ${DATAMASH_DIR}/datamash transpose \
		>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}.QC_REPORT_TEMP.txt
	}

		export -f GRAB_VCF_METRICS_BAIT

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_VCF_METRICS_BAIT \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

########################################################################
##### VCF METRICS FOR TARGET BED FILE ##################################
########################################################################
##### THIS IS THE HEADER ################################################################################
##### "TARGET_HET_HOMVAR_RATIO","TARGET_PCT_GQ0_VARIANT","TARGET_TOTAL_GQ0_VARIANT" #####################
##### "TARGET_TOTAL_HET_DEPTH_SNV","TARGET_TOTAL_SNV","TARGET_NUM_IN_DBSNP_138_SNV" #####################
##### "TARGET_NOVEL_SNV","TARGET_FILTERED_SNV","TARGET_PCT_DBSNP_138_SNV","TARGET_TOTAL_INDEL" ##########
##### "TARGET_NOVEL_INDEL","TARGET_FILTERED_INDEL","TARGET_PCT_DBSNP_138_INDEL" #########################
##### "TARGET_NUM_IN_DBSNP_138_INDEL","TARGET_DBSNP_138_INS_DEL_RATIO","TARGET_NOVEL_INS_DEL_RATIO" #####
##### "TARGET_TOTAL_MULTIALLELIC_SNP","TARGET_NUM_IN_DBSNP_138_MULTIALLELIC" ############################
##### "TARGET_TOTAL_COMPLEX_INDELS,"TARGET_NUM_IN_DBSNP_138_COMPLEX_INDEL" ##############################
##### "TARGET_SNP_REFERENCE_BIAS","TARGET_NUM_SINGLETONS" ###############################################
#########################################################################################################

echo RETRIEVING ON TARGET VCF METRICS: `date`
echo

	GRAB_VCF_METRICS_TARGET ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		awk 'BEGIN {OFS="\t"} \
				$1=="'${SM_TAG}'" \
				{print $2,$3*100,$4,$5,$6,$7,$8,$9,$10*100,$13,$14,$15,$16*100,\
					$17,$18,$19,$20,$21,$22,$23,$24,$25}' \
			${CORE_PATH}/${PROJECT_MS}/REPORTS/VCF_METRICS/MULTI_SAMPLE/${PREFIX}_TARGET.variant_calling_detail_metrics.txt \
			| ${DATAMASH_DIR}/datamash transpose \
		>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}.QC_REPORT_TEMP.txt
	}

		export -f GRAB_VCF_METRICS_TARGET

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_VCF_METRICS_TARGET \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

#########################################################
##### VCF METRICS FOR TITV BED FILE #####################
#########################################################
##### THIS IS THE HEADER ###################################################################################
##### "CODING_HET_HOMVAR_RATIO","CODING_PCT_GQ0_VARIANT","CODING_TOTAL_GQ0_VARIANT" ########################
##### "CODING_TOTAL_HET_DEPTH_SNV","CODING_TOTAL_SNV","CODING_NUM_IN_DBSNP_129_SNV","CODING_NOVEL_SNV" #####
##### "CODING_FILTERED_SNV","CODING_PCT_DBSNP_129_SNV","CODING_DBSNP_129_TITV" #############################
##### "CODING_NOVEL_TITV","CODING_TOTAL_INDEL","CODING_NOVEL_INDEL","CODING_FILTERED_INDEL" ################
##### "CODING_PCT_DBSNP_129_INDEL","CODING_NUM_IN_DBSNP_129_INDEL","CODING_DBSNP_129_INS_DEL_RATIO" ########
##### "CODING_NOVEL_INS_DEL_RATIO","CODING_TOTAL_MULTIALLELIC_SNV" #########################################
##### "CODING_NUM_IN_DBSNP_129_MULTIALLELIC_SNV","CODING_TOTAL_COMPLEX_INDEL" ##############################
##### "CODING_NUM_IN_DBSNP_129_COMPLEX_INDEL","CODING_SNP_REFERENCE_BIAS","CODING_NUM_SINGLETONS" ##########
############################################################################################################

echo RETRIEVING ON TITV VCF METRICS: `date`
echo

	GRAB_VCF_METRICS_TITV ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		awk 'BEGIN {OFS="\t"} \
				$1=="'${SM_TAG}'" \
				{print $2,$3*100,$4,$5,$6,$7,$8,$9,$10*100,$11,$12,$13,$14,$15,$16*100,\
					$17,$18,$19,$20,$21,$22,$23,$24,$25}' \
			${CORE_PATH}/${PROJECT_MS}/REPORTS/VCF_METRICS/MULTI_SAMPLE/${PREFIX}_TITV.variant_calling_detail_metrics.txt \
			| ${DATAMASH_DIR}/datamash transpose \
		>> ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}.QC_REPORT_TEMP.txt
	}

		export -f GRAB_VCF_METRICS_TITV

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		GRAB_VCF_METRICS_TITV \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX

##########################################################################################################
##### tranpose from rows to list so these can be concatenated together for a project/batch QC report #####
##########################################################################################################

echo TRANSPOSING METRICS FROM ROWS TO COLUMNS: `date`
echo

	TRANSPOSE_METRICS ()
	{
		PROJECT_SAMPLE=$1
		SM_TAG=$2
		CORE_PATH=$3
		DATAMASH_DIR=$4
		PROJECT_MS=$5
		PREFIX=$6

		cat ${CORE_PATH}/${PROJECT_MS}/TEMP/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_TEMP.txt" \
		| ${DATAMASH_DIR}/datamash transpose \
		>| ${CORE_PATH}/${PROJECT_MS}/REPORTS/QC_REPORT_PREP_MS/QC_REPORT_PREP_${PREFIX}/${SM_TAG}".QC_REPORT_PREP.txt"
	}

		export -f TRANSPOSE_METRICS

		awk 1 ${SAMPLE_SHEET} \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS=","} NR>1 {print $1,$8}' \
			| sort \
			| uniq \
		| $PARALLEL_DIR/parallel \
			--no-notice \
			-j 4 \
			--colsep ' ' \
		TRANSPOSE_METRICS \
			{1} \
			{2} \
			$CORE_PATH \
			$DATAMASH_DIR \
			$PROJECT_MS \
			$PREFIX
