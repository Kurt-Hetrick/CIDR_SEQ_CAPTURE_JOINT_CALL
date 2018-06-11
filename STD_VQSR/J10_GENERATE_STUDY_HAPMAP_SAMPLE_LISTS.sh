#$ -S /bin/bash
#$ -q rnd.q,prod.q,test.q,bigdata.q
#$ -cwd
#$ -V
#$ -p -50

CORE_PATH=$1
PROJECT=$2
PREFIX=$3

mkdir -p $CORE_PATH/$PROJECT/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/

egrep -m 1 '^#CHROM' $CORE_PATH'/'$PROJECT'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf' | sed 's/\t/\n/g' | awk 'NR>9 {print $0}' | egrep '^NA|^HG' >| $CORE_PATH/$PROJECT/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/$PREFIX'_hapmap_samples.list'

egrep -m 1 '^#CHROM' $CORE_PATH'/'$PROJECT'/MULTI_SAMPLE/'$PREFIX'.BEDsuperset.VQSR.1KG.ExAC3.REFINED.vcf' | sed 's/\t/\n/g' | awk 'NR>9 {print $0}' | egrep -v '^NA|^HG' >| $CORE_PATH/$PROJECT/MULTI_SAMPLE/VARIANT_SUMMARY_STAT_VCF/$PREFIX'_study_samples.list'