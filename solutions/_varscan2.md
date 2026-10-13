```{.bash}
# SS: 1 = Germline, 2 = Somatic, 3 = LOH
module load mugqic/bcftools/1.15
bcftools query -f '%INFO/SS\n' pairedVariants/HCC1395.varscan2.vcf.gz | \
  sort | uniq -c | \
  awk '{s[$2]=$1} END {print s[1]" Germline\n"s[3]" LOH\n"s[2]" Somatic"}'
```
```
526 were called Germline
112 were called LOH
58 were called Somatic
```