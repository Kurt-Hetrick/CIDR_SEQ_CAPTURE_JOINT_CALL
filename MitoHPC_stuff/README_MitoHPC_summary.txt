#### This is the README for files generated from MitoHPC output
 -- het/homo*
 -- per_samp*
 -- all_variants_00*

## Variant level filters applied
 -- FILTER column: no "strict_strand|strand_bias|base_qual|map_qual|weak_evidence|slippage|position|Homopolymer|clustered|fragment|haplotype"
 -- INFO column: no "INDEL|HP(homopolymer region)"
 -- no Read_depth < 300
 -- no variants with missing mlc/mlcm score
 -- no biallelic het variants (for each sample only one het at each POS)

## Sample level filters to apply (not yet applied, filter indicated through reach column)
 -- contamination: contaminated samples flagged from haplocheck
 -- sus_*: suspicious samples from MitoHPC (recommend to filter with sus_mean and sus_fail)
 -- count_het > 5
 -- qc_assessment: delta_CT / CN QC assessment, pass indicated through 'PASSED QC'

## het and homo file
# All het/homo variants passing all variant level filters
 -- filter_numt: whether the variant is NUMT or not
 -- GENE/COMPLEX: the gene / complex info of the variant
 -- mito_lc_consequence/mlc: mito_lc score and consequence for the variant
 -- mlcm: mlc score adjusted with # of homoplasmy seen for each site in UKB (mito_lc_score / (1+log(count_homo+1,10)))
 -- AP_score: Apogee score
 -- MCC_score: missense_OEUF
 -- Median_AF_Het: The median AF of each unique het variant across all samples
 -- Max_AF_Het: The max AF of each unique het variant across all samples
 -- mutation/substitution: nucleotide change and Ti/Tv
 -- mutation_nonsynonymous: whether the mutation is nonsynonymous or not
 -- mutation_stop: nonsense mutation or not
 -- mito_mutation_likelihood: mutation likelihood at each mito location
 -- hom_mlc_score/consequence: homoplasmy version of the mlc score
 -- count_het/homo: number of times each variant is observed in current cohort

## per_sample file
# All samples from MitoHPC pipeline, whether with het or not
 -- count_het/homo: het/homo variants count for each sample
 -- MSS/mMSS: mlc/mlcm sum of all het variants for each sample
 -- contamination: whether the sample is contaminated or not according to haplocheck
 -- CN/delta_ct: mtDNA copy number calculated
 -- haplogroup: haplogroup of each sample
 -- sus_*: suspicious flags from MitoHPC
 -- qc_assessment: CN qc from ALB and DLP
 -- mean: mean mito coverage across mito genome
