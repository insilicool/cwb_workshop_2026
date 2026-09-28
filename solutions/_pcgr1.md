Hunt for another HIGH or MODERATE impact variant.

```{.bash}
less -S pairedVariants/pcgr/HCC1395.pcgr.grch38.pass.vcf.gz
```

There is a MODERATE impact variant in the TP53 (17:7675088) gene, which is a well-known tumor suppressor gene. The variant is a missense mutation that changes an amino acid in the protein sequence.

```{.bash}
# bcftools (use the VEP plugin to extract the CSQ field)
 module purge && \
 module load mugqic/bcftools/1.23 && \
 bcftools +split-vep -r "17:7675088" \
 pairedVariants/pcgr/HCC1395.pcgr.grch38.pass.vcf.gz \
 -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Feature\t%Consequence\t%HGVSc\t%HGVSp\t%CANONICAL\t%MANE_SELECT\n' \
 -d \
 -A tab \
 2>/dev/null
```
There are 38 transcripts for this variant, but the canonical transcript is ENST00000269305.9, which is also the MANE Select transcript (NM_000546.6) and CANONICAL=YES. The consequence of this variant on this transcript is a missense mutation (p.Arg273His), which is classified as a MODERATE impact variant.

```
17	7675088	C	T	TP53	ENST00000269305	missense_variant	ENST00000269305.9:c.524G>A	ENSP00000269305.4:p.Arg175His	YES	NM_000546.6
```

```{.bash}
# bcftools query
module purge && \
 module load mugqic/bcftools/1.23 && \
bcftools query -r "17:7675088" \
-f '%CHROM\t%POS\t%INFO/TSG\t%INFO/TSG_SUPPORT\t%INFO/BIOMARKER_MATCH\t%INFO/ONCOGENICITY\t%INFO/ONCOGENICITY_CODE\n' \
pairedVariants/pcgr/HCC1395.pcgr.grch38.pass.vcf.gz
```
`TSG;TSG_SUPPORT=NCG&CancerMine:1690` - The bare `TSG` flag confirms TP53 is classified as a tumor suppressor gene, sourced from the Network of Cancer Genes (NCG) database and confirmed by CancerMine
The 1690 is CancerMine's citation count — the number of publications supporting TP53's tumor-suppressor role
`ONCOGENICITY=Likely_Oncogenic;ONCOGENICITY_CODE=ONCG_OS1|ONCG_OP1`
OS1 ("Oncogenic Strong-1") — This is the code that fires specifically because p.Arg175His has already been established as oncogenic in prior curated resources.
PCGR isn't evaluating this mutation in isolation, it's recognizing that this exact amino acid substitution has already been classified as oncogenic elsewhere (unlike a novel, never-before-seen missense change, 
which would need to lean on weaker evidence like protein-domain location or in silico prediction alone).
OP1 ("Oncogenic Supporting-1") — "All utilized lines of computational evidence support an oncogenic effect of the variant" (conservation, in silico deleteriousness predictors, etc.). 
This is the same category of dbNSFP-style prediction scores (SIFT, PolyPhen, REVEL, CADD, etc.), here collectively agreeing this missense change is damaging.