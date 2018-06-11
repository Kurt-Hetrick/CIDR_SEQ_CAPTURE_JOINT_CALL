#$ -S /bin/bash
#$ -q rnd.q,prod.q,test.q
#$ -cwd
#$ -p -1000
#$ -V


JAVA_1_7=$1
GATK_DIR=$2
KEY=$3
REF_GENOME=$4

CORE_PATH=$5
PROJECT=$6
SM_TAG=$7

START_SAMPLE_PASS_BAIT_INDEL=`date '+%s'`

CMD=$JAVA_1_7'/java -jar'
CMD=$CMD' '$GATK_DIR'/GenomeAnalysisTK.jar'
CMD=$CMD' -T SelectVariants'
CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'
CMD=$CMD' -et NO_ET'
CMD=$CMD' -K '$KEY
CMD=$CMD' -R '$REF_GENOME
CMD=$CMD' -sn '$SM_TAG
CMD=$CMD' -ef'
CMD=$CMD' -env'
CMD=$CMD' --keepOriginalAC'
CMD=$CMD' -selectType INDEL'
CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT'/VCF/RELEASE/FILTERED_ON_BAIT/'$SM_TAG'_MS_OnBait.vcf'
CMD=$CMD' -o '$CORE_PATH'/'$PROJECT'/INDEL/RELEASE/FILTERED_ON_BAIT/'$SM_TAG'_MS_OnBait_INDEL.vcf'

END_SAMPLE_PASS_BAIT_INDEL=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",L01,SAMPLE_PASS_BAIT_INDEL,"$HOSTNAME","$START_SAMPLE_PASS_BAIT_INDEL","$END_SAMPLE_PASS_BAIT_INDEL \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $CMD >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD | bash
