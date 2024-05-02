# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge that you are in the users current working directory
#$ -cwd

# tell sge to export the users environment variables
#$ -V

# tell sge to submit at this priority setting
#$ -p -10

# tell sge to output both stderr and stdout to the same file
#$ -j y

# export all variables, useful to find out what compute node the program was executed on

	set

	echo

# INPUT VARIABLES

	CORE_PATH=$1

	PROJECT_MS=$2
	SM_TAG=$3
	BAIT_BED=$4
		BAIT_BED_NAME=$(basename ${BAIT_BED} .bed)
	TARGET_BED=$5
		TARGET_BED_NAME=$(basename ${TARGET_BED} .bed)
	TITV_BED=$6
		TITV_BED_NAME=$(basename ${TITV_BED} .bed)
	REF_DICT=$7
	B37_TO_HG19_CHAIN=$8
	HG19_TO_GRCH38_CHAIN=$9


# FIX BED FILES (FOR GRCH37)

	# FIX THE BAIT BED FILE

		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# remove CHR PREFIXES (THIS IS FOR GRCH37)
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
					
			awk 1 ${BAIT_BED} \
				| sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
			>| ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${BAIT_BED_NAME}.bed

	# FIX THE TARGET BED FILE

		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# remove CHR PREFIXES (THIS IS FOR GRCH37)
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
					
			awk 1 ${TARGET_BED} \
				| sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
			>| ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}.bed

	# FIX THE TITV BED FILE

		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# remove CHR PREFIXES (THIS IS FOR GRCH37)
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
					
			awk 1 ${TITV_BED} \
				| sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
			>| ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TITV_BED_NAME}.bed

# # MAKE PICARD INTERVAL FILE (1-based start) for the positions found in the genotyping array vcf if present.

# 	# look for genotyping array vcf and store it as a variable.

# 		GENOTYPING_ARRAY_VCF_TEST=$(ls -tr ${CORE_PATH}/${PROJECT_MS}/Pretesting/Final_Genotyping_VCFs/${SM_TAG}* \
# 			| tail -n 1 )

# # if final report exists containing the full sm-tag, then you're ready to do some formatting.

# if [[ ! -z "${GENOTYPING_ARRAY_VCF_TEST}" ]];then
	
# 		GENOTYPING_ARRAY_VCF=${GENOTYPING_ARRAY_VCF_TEST}

# # if it does not exist, and if the $SM_TAG does not begin with an integer then split $SM_TAG On a @ or -\
# # look for a final report that contains that that first element of the $SM_TAG
# # bonus feature. if this first tests true but the file still does not exist then cidrseqsuite magic files b/c no file exists

# 	elif [[ ${SM_TAG} != [0-9]* ]]; then
		
# 		HAPMAP=${SM_TAG%[@-_]*}
	
# 		GENOTYPING_ARRAY_VCF=$(ls ${CORE_PATH}/${PROJECT_MS}/Pretesting/Final_Genotyping_VCFs/${HAPMAP}* \
# 			| head -n 1)

# else

# # both conditions fails then echo the below message and give a dummy value for the $FINAL_REPORT

# 	echo
# 	echo At this time, you are looking for a genotyping array vcf that does not exist or fails to meet the current logic for finding it.
# 	echo Please talk to Kurt, because he loves to talk.
# 	echo

# 	FINAL_REPORT="FILE_DOES_NOT_EXIST"

# fi

# MAKE PICARD INTERVAL FILES (1-based start) for bed files in the sample sheet
	# GRAB THE SEQUENCING DICTIONARY FORM THE ".dict" file in the directory where the reference genome is located
	# then concatenate with the fixed bed file.
	# add 1 to the start
	# picard interval needs strand information and a locus name
		# made everything plus stranded b/c i don't think this information is used
		# constructed locus name with chr name, start+1, stop

	# bait bed

		(grep "^@SQ" ${REF_DICT} \
			; awk 'BEGIN {OFS="\t"} \
				{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' \
			${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${BAIT_BED_NAME}.bed) \
		>| ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${BAIT_BED_NAME}-picard.bed

	# target bed

		(grep "^@SQ" ${REF_DICT} \
			; awk 'BEGIN {OFS="\t"} \
				{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' \
			${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}.bed) \
		>| ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}-picard.bed

	# titv bed

		(grep "^@SQ" ${REF_DICT} \
			; awk 'BEGIN {OFS="\t"} \
				{print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' \
			${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TITV_BED_NAME}.bed) \
		>| ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TITV_BED_NAME}-picard.bed

# LIFTOVER GRCH37 BED FILE TO HG19 AND THEN GRCH38 TO DO CONCORDANCE TO ARRAY GENOTYPES IF THEY ARE ON GRCH38

	# LIFTOVER TARGET BED FILE TO HG19

		CMD="singularity exec ${ALIGNMENT_CONTAINER} liftOver"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}.bed"
			CMD=${CMD}" ${B37_TO_HG19_CHAIN}"
		CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}-LIFT_HG19.bed"
		CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}-LIFT_HG19_REJECTED.bed"

		# write command line to file and execute the command line

			echo ${CMD} >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${SM_TAG}_command_lines.txt
			echo >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${SM_TAG}_command_lines.txt
			echo ${CMD} | bash

	# LIFTOVER HG19 BED FILE TO GRCH38

		CMD="singularity exec ${ALIGNMENT_CONTAINER} liftOver"
			CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}-LIFT_HG19.bed"
			CMD=${CMD}" ${HG19_TO_GRCH38_CHAIN}"
		CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}-LIFT_GRCH38.bed"
		CMD=${CMD}" ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}-LIFT_GRCH38_REJECTED.bed"

		# write command line to file and execute the command line

			echo ${CMD} >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${SM_TAG}_command_lines.txt
			echo >> ${CORE_PATH}/${PROJECT_MS}/COMMAND_LINES/${SM_TAG}_command_lines.txt
			echo ${CMD} | bash

	# remove any loci that are not part of the primary assembly
	# this is for concordance when the gt array reference genome is grch38 b/c cidrseqsuite will crash

		grep -v "^@" ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}-LIFT_GRCH38.bed \
			| awk 'BEGIN {OFS="\t"} \
				$1!~"_" \
				{print $0}' \
		>| ${CORE_PATH}/${PROJECT_MS}/TEMP/${SM_TAG}/${SM_TAG}-${TARGET_BED_NAME}-LIFT_GRCH38_PRIMARY.bed
