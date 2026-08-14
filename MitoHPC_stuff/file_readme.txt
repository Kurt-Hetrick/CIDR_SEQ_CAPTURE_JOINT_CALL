count.tab: total reads and mtDNA-CN counts
cvg.tab: coverage stats
het.txt : All het variants passing all variant level filters. README_MitoHPC_summary.txt
homo.txt : All homo variants passing all variant level filters. README_MitoHPC_summary.txt
per_sample.txt: sample summary stats. see README_MitoHPC_summary.txt

mutect2.cvg.tab : coverage stats using sample's consensus sequence
mutect2.fa : consensus sequence for sample's mitochondrial genome
mutect2.fa.fai : index file for mutect2.fa
mutect2.haplocheck.tab : contamination screen   
mutect2.haplogroup1.tab : major haplogroup
mutect2.haplogroup.tab : haplogroup
mutect2.merge.bed

mutect2.00.AF.histo : histrogram of allele frequency
mutect2.00.concat.vcf : single SNV from single sample in each line; GT/DP/AF

=======================
filtered to > 3% VAF
=======================
mutect2.03.vcf : SNVs; 3% heteroplasmy thold; single SNV from single sample in each line; last column = sample name
mutect2.03.concat.vcf : SNVs; 3% heteroplasmy thold; single SNV from single sample in each line; GT/DP/AF
mutect2.03.merge.vcf : SNVs; 3% heteroplasmy threshold; single SNV from multiple samples in each line
mutect2.03.merge.sitesOnly.vcf : SNVs; 3% heteroplasmy threshold; single SNV from multiple samples in each line, excluding sample specific info
mutect2.03.pos : position summaries; AC and AF values for both HOM and HET
mutect2.03.summary : SNV count summaries
mutect2.03.suspicious.ids : samples which either have low mtDNA-CN,low cvg, failed haplockeck, multiple NUMT/HG labels (AF<1)
mutect2.03.suspicious.tab : samples which either have low mtDNA-CN,low cvg, failed haplockeck, multiple NUMT/HG labels (AF<1) (one line per sample/reason, can have multiple lines for a sample)
mutect2.03.tab : SNV counts

=======================
filtered to > 5% VAF
=======================
mutect2.05.vcf : SNVs; 5% heteroplasmy thold; single SNV from single sample in each line; last column = sample name
mutect2.05.concat.vcf : SNVs; 5% heteroplasmy thold; single SNV from single sample in	each line; GT/DP/AF
mutect2.05.merge.vcf : SNVs; 5% heteroplasmy threshold; single SNV from multiple samples in each line
mutect2.05.merge.sitesOnly.vcf : SNVs; 5% heteroplasmy threshold; single SNV from multiple samples in each line, excluding sample specific info
mutect2.05.pos : position summaries; AC and AF values for both HOM and HET
mutect2.05.summary : SNV count summaries
mutect2.05.suspicious.ids : samples which either have low mtDNA-CN,low cvg, failed haplockeck, multiple NUMT/HG labels (AF<1)
mutect2.05.suspicious.tab : samples which either have low mtDNA-CN,low cvg, failed haplockeck, multiple NUMT/HG labels (AF<1) (one line per sample/reason, can have multiple lines for a sample)
mutect2.05.tab : SNV counts

=======================
filtered to > 10% VAF
=======================
mutect2.10.vcf : SNVs; 10% heteroplasmy thold; single SNV from single sample in each line; last column = sample name
mutect2.10.concat.vcf : SNVs; 10% heteroplasmy thold; single SNV from single sample in	each line; GT/DP/AF
mutect2.10.merge.vcf : SNVs; 10% heteroplasmy threshold; single SNV from multiple samples in each line
mutect2.10.merge.sitesOnly.vcf : SNVs; 10% heteroplasmy threshold; single SNV from multiple samples in each line, excluding sample specific info
mutect2.10.pos : position summaries; AC and AF values for both HOM and HET
mutect2.10.summary : SNV count summaries
mutect2.10.suspicious.ids : samples which either have low mtDNA-CN,low cvg, failed haplockeck, multiple NUMT/HG labels (AF<1)
mutect2.10.suspicious.tab : samples which either have low mtDNA-CN,low cvg, failed haplockeck, multiple NUMT/HG labels (AF<1) (one line per sample/reason, can have multiple lines for a sample)
mutect2.10.tab : SNV counts
===============================



=======================
these are the 2nd iteration calls. variants are called against a sample's consensus sequence.
=======================
mutect2.mutect2.00.AF.histo : histrogram of allele frequency
mutect2.mutect2.00.concat.vcf : single SNV from single sample in each line; GT/DP/AF

=======================
filtered to > 3% VAF
=======================
mutect2.mutect2.03.vcf : SNVs; 3% heteroplasmy thold; single SNV from single sample in each line; last column = sample name
mutect2.mutect2.03.concat.vcf : SNVs; 3% heteroplasmy thold; single SNV from single sample in each line; GT/DP/AF
mutect2.mutect2.03.merge.vcf : SNVs; 3% heteroplasmy threshold; single SNV from multiple samples in each line
mutect2.mutect2.03.merge.sitesOnly.vcf : SNVs; 3% heteroplasmy threshold; single SNV from multiple samples in each line, excluding sample specific info
mutect2.mutect2.03.pos : position summaries; AC and AF values
mutect2.mutect2.03.summary : SNV count summaries
mutect2.mutect2.03.tab : SNV counts

=======================
filtered to > 5% VAF
=======================
mutect2.mutect2.05.vcf : SNVs; 5% heteroplasmy thold; single SNV from single sample in each line; last column = sample name
mutect2.mutect2.05.concat.vcf : SNVs; 5% heteroplasmy thold; single SNV from single sample in each line; GT/DP/AF
mutect2.mutect2.05.merge.vcf : SNVs; 5% heteroplasmy threshold; single SNV from multiple samples in each line
mutect2.mutect2.05.merge.sitesOnly.vcf : SNVs; 5% heteroplasmy threshold; single SNV from multiple samples in each line, excluding sample specific info
mutect2.mutect2.05.pos : position summaries; AC and AF values
mutect2.mutect2.05.summary : SNV count summaries
mutect2.mutect2.05.tab : SNV counts

=======================
filtered to > 10% VAF
=======================
mutect2.mutect2.10.vcf : SNVs; 10% heteroplasmy thold; single SNV from single sample in each line; last column = sample name
mutect2.mutect2.10.concat.vcf : SNVs; 10% heteroplasmy thold; single SNV from single sample in each line; GT/DP/AF
mutect2.mutect2.10.merge.vcf : SNVs; 10% heteroplasmy threshold; single SNV from multiple samples in each line
mutect2.mutect2.10.merge.sitesOnly.vcf : SNVs; 10% heteroplasmy threshold; single SNV from multiple samples in each line, excluding sample specific info
mutect2.mutect2.10.pos : position summaries; AC and AF values
mutect2.mutect2.10.summary : SNV count summaries
mutect2.mutect2.10.tab : SNV counts
