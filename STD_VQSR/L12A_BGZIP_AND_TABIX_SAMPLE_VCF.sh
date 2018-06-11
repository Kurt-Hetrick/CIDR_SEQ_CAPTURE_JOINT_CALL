#$ -S /bin/bash
#$ -q rnd.q,prod.q,test.q
#$ -cwd
#$ -p -1000
#$ -V

TABIX_DIR=$1

CORE_PATH=$2
PROJECT=$3
SM_TAG=$4

START_BGZIP_AND_TABIX=`date '+%s'`

CMD1=$TABIX_DIR'/bgzip -c '$CORE_PATH'/'$PROJECT'/VCF/RELEASE/FILTERED_ON_BAIT/'$SM_TAG'_MS_OnBait.vcf'
CMD1=$CMD1' >| '$CORE_PATH'/'$PROJECT'/VCF/RELEASE/FILTERED_ON_BAIT/TABIX/'$SM_TAG'_MS_OnBait.vcf.gz'

CMD2=$TABIX_DIR'/tabix -f -p vcf '$CORE_PATH'/'$PROJECT'/VCF/RELEASE/FILTERED_ON_BAIT/TABIX/'$SM_TAG'_MS_OnBait.vcf.gz'

END_BGZIP_AND_TABIX=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",L01,BGZIP_AND_TABIX,"$HOSTNAME","$START_BGZIP_AND_TABIX","$END_BGZIP_AND_TABIX \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $CMD1 >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD1 | bash

echo $CMD2 >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD2 | bash

