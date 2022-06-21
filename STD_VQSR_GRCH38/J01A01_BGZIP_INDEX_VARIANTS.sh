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

	TABIX_DIR=$1

	CORE_PATH=$2
	PROJECT_MS=$3
	PREFIX=$4

START_BGZIP_INDEX=`date '+%s'`

# compress vcf file with bgzip

	CMD1=$TABIX_DIR'/bgzip -c '$CORE_PATH'/'$PROJECT_MS'/TEMP/'$PREFIX'.BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf'
	CMD1=$CMD1' >| '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf.gz'

# index vcf file

	CMD2=$TABIX_DIR'/tabix -p vcf -f '$CORE_PATH'/'$PROJECT_MS'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf.gz'

END_BGZIP_INDEX=`date '+%s'`

# send command line to command line file

	# bgzip command

		echo $CMD1 >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
		echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
		echo $CMD1 | bash

	# tabix command

		echo $CMD2 >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
		echo >> $CORE_PATH/$PROJECT_MS/COMMAND_LINES/$PROJECT_MS"_command_lines.txt"
		echo $CMD2 | bash

echo $PROJECT_MS",K01,BGZIP_VARIANTS,"$HOSTNAME","$START_BGZIP_INDEX","$END_BGZIP_INDEX \
>> $CORE_PATH/$PROJECT_MS/REPORTS/$PROJECT_MS".JOINT.CALL.WALL.CLOCK.TIMES.csv"

# check to see if the index is generated which should send an non-zero exit signal if not.
# eventually, will want to check the exit signal above and push out whatever it is at the end. Not doing that today though.

ls $CORE_PATH/$PROJECT_MS/MULTI_SAMPLE/$PREFIX".BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf.gz.tbi"
