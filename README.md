# PopPsiSeq
Population-based framework for introgression/selection/resequencing experiments

# Background
Designed to be an extension to the "PsiSeq" protocol, which identifies candidate genome locii in selection/backcross experiments (Earley & Jones 2011).


# What's New
PopPsiSeq updates the original protocol in several ways:

* consolidation into a Snakemake pipeline simplifies a somewhat unruly workflow
* improved data QA/QC, quality & uniqueness of mapping
* uses empirical sequenced reads rather than fragmented reference genome to characterize species.
* incorporates software advances in variant calling (eg, freebayes rather than directly examining pileup), data processing & visualization (eg, ggplot), and other utilities (eg, vcftools, bedtools) 
* PsiSeq uses a reciprocal mapping scheme to call variants (eg, simulans reads vs sechellia reference and sechellia reads vs simulans reference), whereas PopPsiSeq currently maps both to a third, common reference (eg, simulans reads &  sechellia reads vs melanogaster reference).
* PsiSeq assumes that differences between species are fixed; PopPsiSeq examines local changes in allele frequency (of which fixation is an extreme case). 

# References

Earley, Eric J., and Corbin D. Jones. 2011. “Next-generation mapping of complex traits with phenotype-based selection and introgression.” Genetics 189 (4): 1203–9. doi:10.1534/genetics.111.129445.



