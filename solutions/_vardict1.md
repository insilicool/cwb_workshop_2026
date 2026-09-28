The other Vardict variant were filtered out due to the -f PASS argument in bcftools view command by looking at the FILTER column in the VCF file. 
The FILTER column indicates whether a variant passed or failed the filtering criteria.

The reasons for filtering can be found in the ##FILTER lines in the VCF header. A variant may fail one or more of these filters, and the specific filter(s) that caused the variant to be filtered out will be listed in the FILTER column for that variant.

```{.bash}
zgrep "##FILTER" pairedVariants/HCC1395.vardict.vcf.gz
```

```
##FILTER=<ID=q22.5,Description="Mean Base Quality Below 22.5">
##FILTER=<ID=Q0,Description="Mean Mapping Quality Below 0">
##FILTER=<ID=p8,Description="Mean Position in Reads Less than 8">
##FILTER=<ID=SN1.5,Description="Signal to Noise Less than 1.5">
##FILTER=<ID=Bias,Description="Strand Bias">
##FILTER=<ID=pSTD,Description="Position in Reads has STD of 0">
##FILTER=<ID=MAF0.05,Description="Matched sample has AF > 0.05, thus not somatic">
##FILTER=<ID=d5,Description="Total Depth < 5">
##FILTER=<ID=v3,Description="Var Depth < 3">
##FILTER=<ID=f0.03,Description="Allele frequency < 0.03">
##FILTER=<ID=P0.9,Description="Not significant with p-value > 0.9">
##FILTER=<ID=DIFF0.2,Description="Non-somatic or LOH and allele frequency difference < 0.2">
##FILTER=<ID=P0.01Likely,Description="Likely candidate but p-value > 0.01/5**vd2">
##FILTER=<ID=InDelLikely,Description="Likely Indels are not considered somatic">
##FILTER=<ID=MSI12,Description="Variant in MSI region with 12 non-monomer MSI or 12 monomer MSI">
##FILTER=<ID=NM4.25,Description="Mean mismatches in reads >= 4.25, thus likely false positive">
##FILTER=<ID=InGap,Description="The somatic variant is in the deletion gap, thus likely false positive">
##FILTER=<ID=InIns,Description="The somatic variant is adjacent to an insertion variant">
##FILTER=<ID=Cluster0bp,Description="Two somatic variants are within 0 bp">
##FILTER=<ID=LongAT,Description="The somatic variant is flanked by long A/T (>=14)">

```

Also it should be noted that selecting by STATUS also removes variants denoted as AFDiff and Deletion, which are not considered somatic variants.  Variants with other STATUS values, such as "AFDiff" or "Deletion", are filtered out and not included in the final set of somatic and germline/LOH variants.