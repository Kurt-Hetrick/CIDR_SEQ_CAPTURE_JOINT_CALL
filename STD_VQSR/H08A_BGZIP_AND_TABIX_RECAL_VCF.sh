#$ -S /bin/bash
#$ -q rnd.q,prod.q,test.q
#$ -cwd
#$ -V
#$ -p -1000

TABIX_DIR=$1

CORE_PATH=$2
PROJECT=$3
PREFIX=$4

START_BGZIP_TABIX=`date '+%s'`

CMD1=$TABIX_DIR'/bgzip -c '$CORE_PATH'/'$PROJECT'/MULTI_SAMPLE/'$PREFIX'.HC.SNP.INDEL.VQSR.vcf'
CMD1=$CMD1' >| '$CORE_PATH'/'$PROJECT'/MULTI_SAMPLE/'$PREFIX'.HC.SNP.INDEL.VQSR.vcf.gz'

CMD2=$TABIX_DIR'/tabix -p vcf -f '$CORE_PATH'/'$PROJECT'/MULTI_SAMPLE/'$PREFIX'.HC.SNP.INDEL.VQSR.vcf.gz'

END_BGZIP_TABIX=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",H01,BGZIP_TABIX,"$HOSTNAME","$START_BGZIP_TABIX","$END_BGZIP_TABIX \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $CMD1 >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD1 | bash

echo $CMD2 >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD2 | bash
