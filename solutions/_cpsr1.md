format2pcgr.py added:

```
##INFO=<ID=TAL,Number=1,Type=String,Description="Confidence of call i.e. number of callers supporting the variant">
##INFO=<ID=TDP,Number=1,Type=Integer,Description="Tumor depth derived from tumor AD field">
##INFO=<ID=TVAF,Number=1,Type=Float,Description="Tumor variant allele frequency derived from tumor AD field">
##INFO=<ID=NDP,Number=1,Type=Integer,Description="Normal depth derived from tumor AD field">
##INFO=<ID=NVAF,Number=1,Type=Float,Description="Normal variant allele frequency derived from tumor AD field">
```

bcftools added:

```
##bcftools_viewCommand=view -Oz '-iTDP>=10 && TVAF>=0.05 && NDP>=10 && NVAF>=0.05' pairedVariants/HCC1395.ensemble.germline.format.vcf.gz; Date=Thu Sep 24 14:28:44 2026
```
