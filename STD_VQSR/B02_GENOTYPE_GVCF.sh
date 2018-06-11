#$ -S /bin/bash
#$ -q rnd.q,prod.q,test.q
#$ -cwd
#$ -V
#$ -p -1000

JAVA_1_7=$1
GATK_DIR=$2
REF_GENOME=$3
KEY=$4

CORE_PATH=$5
PROJECT=$6
PREFIX=$7
BED_FILE_NAME=$8

START_GENOTYPE_GVCF=`date '+%s'`

CMD=$JAVA_1_7'/java -jar'
CMD=$CMD' '$GATK_DIR'/GenomeAnalysisTK.jar'
CMD=$CMD' -T GenotypeGVCFs'
CMD=$CMD' -R '$REF_GENOME
CMD=$CMD' --annotateNDA'
CMD=$CMD' --variant '$CORE_PATH'/'$PROJECT'/GVCF/AGGREGATE/'$PREFIX'.'$BED_FILE_NAME'.genome.vcf'
CMD=$CMD' --disable_auto_index_creation_and_locking_when_reading_rods'
CMD=$CMD' -et NO_ET'
CMD=$CMD' -K '$KEY
CMD=$CMD' -o '$CORE_PATH'/'$PROJECT'/TEMP/'$PREFIX'.'$BED_FILE_NAME'.temp.vcf'

echo $CMD >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD | bash

END_GENOTYPE_GVCF=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",B01,GENOTYPE_GVCF,"$HOSTNAME","$START_GENOTYPE_GVCF","$END_GENOTYPE_GVCF \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"
