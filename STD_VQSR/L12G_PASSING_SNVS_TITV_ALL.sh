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
TITV_BED=$8

START_ALL_SNV_TITV=`date '+%s'`

CMD=$JAVA_1_7'/java -jar'
CMD=$CMD' '$GATK_DIR'/GenomeAnalysisTK.jar'
CMD=$CMD' -T SelectVariants'
CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'
CMD=$CMD' -et NO_ET'
CMD=$CMD' -K '$KEY
CMD=$CMD' -R '$REF_GENOME
CMD=$CMD' -ef'
CMD=$CMD' -env'
CMD=$CMD' --keepOriginalAC'
CMD=$CMD' -L '$TITV_BED
CMD=$CMD' -selectType SNP'
CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT'/VCF/RELEASE/FILTERED_ON_BAIT/'$SM_TAG'_MS_OnBait.vcf'
CMD=$CMD' -o '$CORE_PATH'/'$PROJECT'/TEMP/'$SM_TAG'.Release.OnExon.FILTERED.vcf'

END_ALL_SNV_TITV=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",L01,ALL_SNV_TITV,"$HOSTNAME","$START_ALL_SNV_TITV","$END_ALL_SNV_TITV \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $CMD >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD | bash
