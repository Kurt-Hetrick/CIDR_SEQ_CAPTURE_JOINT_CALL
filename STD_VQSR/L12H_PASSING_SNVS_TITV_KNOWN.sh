#$ -S /bin/bash
#$ -q rnd.q,prod.q,test.q
#$ -cwd
#$ -p -1000
#$ -V


JAVA_1_7=$1
GATK_DIR=$2
KEY=$3
REF_GENOME=$4
KNOWN_SNPS=$5

CORE_PATH=$6
PROJECT=$7
SM_TAG=$8
TITV_BED=$9

START_KNOWN_SNV_TITV=`date '+%s'`

CMD=$JAVA_1_7'/java -jar'
CMD=$CMD' '$GATK_DIR'/GenomeAnalysisTK.jar'
CMD=$CMD' -T SelectVariants'
CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'
CMD=$CMD' -et NO_ET'
CMD=$CMD' -K '$KEY
CMD=$CMD' -R '$REF_GENOME
CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT'/VCF/RELEASE/FILTERED_ON_BAIT/'$SM_TAG'_MS_OnBait.vcf'
CMD=$CMD' --excludeFiltered'
CMD=$CMD' --excludeNonVariants'
CMD=$CMD' --keepOriginalAC'
CMD=$CMD' -L '$TITV_BED
CMD=$CMD' --concordance '$KNOWN_SNPS
CMD=$CMD' -selectType SNP'
CMD=$CMD' -o '$CORE_PATH'/'$PROJECT'/TEMP/'$SM_TAG'.Release.Known.TiTv.vcf'

END_KNOWN_SNV_TITV=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",L01,KNOWN_SNV_TITV,"$HOSTNAME","$START_KNOWN_SNV_TITV","$END_KNOWN_SNV_TITV \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $CMD >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD | bash
