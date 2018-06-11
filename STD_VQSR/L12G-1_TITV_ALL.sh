#$ -S /bin/bash
#$ -q rnd.q,prod.q,test.q
#$ -cwd
#$ -p -1000
#$ -V

SAMTOOLS_DIR=$1

CORE_PATH=$2
PROJECT=$3
SM_TAG=$4

START_ALL_TITV=`date '+%s'`

CMD=$SAMTOOLS_DIR'/bcftools/vcfutils.pl qstats '$CORE_PATH'/'$PROJECT'/TEMP/'$SM_TAG'.Release.OnExon.FILTERED.vcf'
CMD=$CMD' >| '$CORE_PATH'/'$PROJECT'/REPORTS/TI_TV_MS/'$SM_TAG'_All_.titv.txt'

END_ALL_TITV=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",M01,ALL_TITV,"$HOSTNAME","$START_ALL_TITV","$END_ALL_TITV \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $CMD >> $CORE_PATH/$PROJECT/command_lines.txt
echo >> $CORE_PATH/$PROJECT/command_lines.txt
echo $CMD | bash
