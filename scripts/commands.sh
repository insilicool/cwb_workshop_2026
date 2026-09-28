#!/usr/bin/env bash
#
# scripts/commands.sh
#
# Auto-derived from README.md: this file concatenates every ```{.bash}``` code
# block in README.md, in document order, into one non-interactive script so
# the whole practical can be run end-to-end without manual intervention.
#
# Regenerate this file whenever README.md's command blocks change.
#
# Assumptions:
#   - You are already running inside the workshop environment (either the
#     c3genomics/genpipes docker container or an HPC session with the
#     mugqic module system available). The `docker run` command that starts
#     that environment is intentionally NOT executed by this script (see
#     block 1 below) -- start the container/session yourself, then run this
#     script from inside it.
#   - $REF, $COURSE and friends are exported by block 2 below; downstream
#     tool paths (e.g. $VARSCAN2_JAR, $VARDICT_HOME, $PCGR_DATA, ...) are
#     provided by the "module load" commands via the mugqic module system.
#
set -euo pipefail

# ===== Block 1/32 =====
# NOTE: run this manually to enter the container/environment BEFORE
# running the rest of this script. Not executed here.
# #Spin up docker container
# docker run --privileged -v /tmp:/tmp --network host -it \
#     -w $PWD -v $HOME:$HOME -v /etc/fonts/:/etc/fonts/ \
#     -v $HOME/cvmfs_caches/:/cvmfs-cache/ c3genomics/genpipes:v6.2.0
#
# #Load modules and reference resources    
# module purge
#
# export REF=${MUGQIC_INSTALL_HOME}/genomes/species/Homo_sapiens.GRCh38/
# export COURSE=/home/training/cbw_workshop_2026
#
# #create and move to working directory
# mkdir -p $COURSE/SNV
#
# cd $COURSE/SNV

# ===== Block 2/32 =====
ls -l alignment/HCC1395BL_normal/
samtools view -H alignment/HCC1395BL_normal/HCC1395BL_normal.subset.bam

# ===== Block 3/32 =====
samtools view -H alignment/HCC1395BL_normal/HCC1395BL_normal.subset.bam | grep "@RG"

# ===== Block 4/32 =====
mkdir -p pairedVariants

# ===== Block 5/32 =====
# SAMTools mpileup
##Good practice to purge all modules before loading new ones
module purge && \
module load mugqic/samtools/1.14 && \
samtools mpileup -d 1000 -B -q 10 -Q 10 \
  -f ${REF}/genome/Homo_sapiens.GRCh38.fa \
  -l regions.bed \
  alignment/HCC1395BL_normal/HCC1395BL_normal.subset.bam \
  alignment/HCC1395_tumor/HCC1395_tumor.subset.bam \
  > pairedVariants/HCC1395.varscan2.mpileup

# ===== Block 6/32 =====
# Varscan2
module purge && \
module load mugqic/java/openjdk-jdk-17.0.1 mugqic/VarScan/2.4.3 mugqic/htslib/1.14 mugqic/bcftools/1.15 && \
java -Xmx2000M -jar $VARSCAN2_JAR somatic \
  pairedVariants/HCC1395.varscan2.mpileup \
  pairedVariants/HCC1395.varscan2 \
  --min-coverage 3 \
  --min-var-freq 0.05 \
  --p-value 0.10 \
  --somatic-p-value 0.05 \
  --strand-filter 0 \
  --output-vcf 1 \
  --mpileup 1 && \
bgzip -cf  \
 pairedVariants/HCC1395.varscan2.snp.vcf \
  > pairedVariants/HCC1395.varscan2.snp.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.varscan2.snp.vcf.gz && \
bgzip -cf  \
 pairedVariants/HCC1395.varscan2.indel.vcf \
  > pairedVariants/HCC1395.varscan2.indel.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.varscan2.indel.vcf.gz && \
bcftools \
  concat -a \
  pairedVariants/HCC1395.varscan2.snp.vcf.gz \
  pairedVariants/HCC1395.varscan2.indel.vcf.gz | \
sed 's/TUMOR/HCC1395_NS_T_1/g'   | \
sed 's/NORMAL/HCC1395BL_NS_N_1/g' \
> pairedVariants/HCC1395.varscan2.vcf && \
bgzip -cf  \
 pairedVariants/HCC1395.varscan2.vcf \
  > pairedVariants/HCC1395.varscan2.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.varscan2.vcf.gz

# ===== Block 7/32 =====
# View Mandatory fields (-E means to use extended regex)
zgrep -E "##fileformat|#CHROM" pairedVariants/HCC1395.varscan2.vcf.gz

# ===== Block 8/32 =====
# View FORMAT fields
zgrep -E "##FORMAT" pairedVariants/HCC1395.varscan2.vcf.gz

# ===== Block 9/32 =====
# View INFO fields
zgrep "##INFO" pairedVariants/HCC1395.varscan2.vcf.gz

# ===== Block 10/32 =====
# Extract somatic variants
module purge && \
module load mugqic/bcftools/1.15 mugqic/htslib/1.14 && \
bcftools view -Oz -i 'INFO/SS="2"' -o pairedVariants/varscan2.somatic.vcf.gz \
pairedVariants/HCC1395.varscan2.vcf.gz && tabix -pvcf pairedVariants/varscan2.somatic.vcf.gz

# Sanity check (remove the header and count the number of lines)
zgrep -v "^#" pairedVariants/varscan2.somatic.vcf.gz | wc -l

# ===== Block 11/32 =====
# View the first somatic variant (-A1 means to show 1 line after the match)
zgrep -A1 "#CHROM" pairedVariants/HCC1395.varscan2.vcf.gz

# ===== Block 12/32 =====
# Extract germline and LOH variants
module purge && \
module load mugqic/bcftools/1.15 mugqic/htslib/1.14 && \
bcftools view -Oz -i 'INFO/SS="1" | INFO/SS="3"' -o pairedVariants/HCC1395.varscan2.germline.loh.vcf.gz \
pairedVariants/HCC1395.varscan2.vcf.gz && tabix -pvcf pairedVariants/HCC1395.varscan2.germline.loh.vcf.gz

##Sanity check (remove the header and count the number of lines)
zgrep -v "^#" pairedVariants/HCC1395.varscan2.germline.loh.vcf.gz | wc -l 

# ===== Block 13/32 =====
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/VarDictJava/1.4.8 mugqic/samtools/1.14 mugqic/perl/5.34.0 mugqic/R_Bioconductor/4.1.0_3.13 mugqic/htslib/1.14 && \
java -Xmx6000M -classpath $VARDICT_HOME/lib/VarDict-1.4.8.jar:$VARDICT_HOME/lib/commons-cli-1.2.jar:$VARDICT_HOME/lib/jregex-1.2_01.jar:$VARDICT_HOME/lib/htsjdk-2.8.0.jar com.astrazeneca.vardict.Main \
  --G ${REF}/genome/Homo_sapiens.GRCh38.fa \
  -N HCC1395 \
  -b "alignment/HCC1395_tumor/HCC1395_tumor.subset.bam|alignment/HCC1395BL_normal/HCC1395BL_normal.subset.bam" \
  -f 0.03 -Q 10 -c 1 -S 2 -E 3 -g 4 -th 3 \
  regions.bed | \
$VARDICT_BIN/testsomatic.R  | \
perl $VARDICT_BIN/var2vcf_paired.pl \
    -N "HCC1395_NS_T_1|HCC1395BL_NS_N_1" \
    -f 0.03 -P 0.9 -m 4.25 -M | \
bgzip -cf  \
  > pairedVariants/HCC1395.vardict.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.vardict.vcf.gz && \
zgrep -v "^#" pairedVariants/HCC1395.vardict.vcf.gz | wc -l

# ===== Block 14/32 =====
## Extract somatic variants
module purge && \
module load mugqic/bcftools/1.15 mugqic/htslib/1.14 && \
bcftools \
  view -f PASS -i 'INFO/STATUS~".*Somatic"' \
  pairedVariants/HCC1395.vardict.vcf.gz | \
bgzip -cf  \
  > pairedVariants/HCC1395.vardict.somatic.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.vardict.somatic.vcf.gz && \
zgrep -v "^#" pairedVariants/HCC1395.vardict.somatic.vcf.gz | wc -l

## Extract germline and LOH variants
module purge && \
module load mugqic/bcftools/1.15  mugqic/htslib/1.14 && \
bcftools \
  view -f PASS -i 'INFO/STATUS~"Germline" | INFO/STATUS~".*LOH"' \
  pairedVariants/HCC1395.vardict.vcf.gz | \
bgzip -cf  \
  > pairedVariants/HCC1395.vardict.germline.loh.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.vardict.germline.loh.vcf.gz && \
zgrep -v "^#" pairedVariants/HCC1395.vardict.germline.loh.vcf.gz | wc -l

# ===== Block 15/32 =====
# Variants MuTecT2
module purge && \
module load mugqic/java/openjdk-jdk-17.0.1 mugqic/GenomeAnalysisTK/4.6.0.0 && \
gatk --java-options "-Xmx6000M" \
  Mutect2 \
  --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP --native-pair-hmm-threads 3 \
  --max-reads-per-alignment-start 0 --read-validation-stringency LENIENT \
  --af-of-alleles-not-in-resource 0.0000025 \
  --f1r2-tar-gz pairedVariants/HCC1395.f1r2.tar.gz \
  --reference ${REF}/genome/Homo_sapiens.GRCh38.fa \
  --input alignment/HCC1395_tumor/HCC1395_tumor.subset.bam \
  --tumor-sample HCC1395_NS_T_1 \
  --input alignment/HCC1395BL_normal/HCC1395BL_normal.subset.bam \
  --normal-sample HCC1395BL_NS_N_1 \
  --germline-resource Homo_sapiens.GRCh38.af-only-gnomad.regions_only.vcf.gz \
  --intervals regions.bed \
  --output pairedVariants/HCC1395.mutect2.vcf.gz

# ===== Block 16/32 =====
# Filtering
module purge && \
module load mugqic/java/openjdk-jdk-17.0.1 mugqic/GenomeAnalysisTK/4.6.0.0 mugqic/bcftools/1.15 mugqic/htslib/1.14 && \
gatk --java-options "-Xmx1000M" \
  LearnReadOrientationModel  \
  --input pairedVariants/HCC1395.f1r2.tar.gz \
  --output pairedVariants/HCC1395.read-orientation-model.tar.gz && \
gatk --java-options "-Xmx1000M" \
  FilterMutectCalls  \
  --reference ${REF}/genome/Homo_sapiens.GRCh38.fa \
  --variant pairedVariants/HCC1395.mutect2.vcf.gz \
  --ob-priors pairedVariants/HCC1395.read-orientation-model.tar.gz \
  --output pairedVariants/HCC1395.mutect2.flt.vcf.gz && \
bcftools view -f PASS -Oz -o pairedVariants/HCC1395.mutect2.somatic.vcf.gz pairedVariants/HCC1395.mutect2.flt.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.mutect2.somatic.vcf.gz && \
zgrep -v "^#" pairedVariants/HCC1395.mutect2.somatic.vcf.gz | wc -l

# ===== Block 17/32 =====
# Strelka2 somatic (region.bed must be bgzipped and tabix indexed)
module purge && \
module load mugqic/htslib/1.14 mugqic/python/2.7.18 mugqic/Strelka2/2.9.10 && \
cat regions.bed | bgzip -cf > regions.bed.gz && \
tabix -f -pbed regions.bed.gz && \
rm -rf pairedVariants/strelka2_somatic && \
python $STRELKA2_HOME/bin/configureStrelkaSomaticWorkflow.py \
  --normalBam alignment/HCC1395BL_normal/HCC1395BL_normal.subset.bam \
  --tumorBam alignment/HCC1395_tumor/HCC1395_tumor.subset.bam \
  --referenceFasta ${REF}/genome/Homo_sapiens.GRCh38.fa \
  --callRegions regions.bed.gz \
  --runDir pairedVariants/strelka2_somatic && \
python pairedVariants/strelka2_somatic/runWorkflow.py \
  -m local  \
  -j 3 \
  -g 6 \
  --quiet

# ===== Block 18/32 =====
# Filter for somatic (update_genotypes_strelka.py adds FORMAT GT to samples; strelka2 uses a non-standard format)
module purge && \
module load mugqic/bcftools/1.15 mugqic/htslib/1.14 mugqic/mugqic_tools/2.12.7 mugqic/python/3.10.4 && \
bcftools \
  concat -a  \
  pairedVariants/strelka2_somatic/results/variants/somatic.snvs.vcf.gz \
  pairedVariants/strelka2_somatic/results/variants/somatic.indels.vcf.gz | \
sed 's/TUMOR/HCC1395_NS_T_1/g'   | \
sed 's/NORMAL/HCC1395BL_NS_N_1/g'   | \
bgzip -cf  \
  > pairedVariants/HCC1395.strelka2.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.strelka2.vcf.gz && \
  python3 $PYTHON_TOOLS/update_genotypes_strelka.py \
      -i pairedVariants/HCC1395.strelka2.vcf.gz \
      -o pairedVariants/HCC1395.strelka2.gt.vcf.gz \
      -n HCC1395BL_NS_N_1 \
      -t HCC1395_NS_T_1 && \
bcftools \
  view -f PASS -Oz \
 -o pairedVariants/HCC1395.strelka2.somatic.vcf.gz \
 pairedVariants/HCC1395.strelka2.gt.vcf.gz

# ===== Block 19/32 =====
module purge && \
module load mugqic/htslib/1.14 mugqic/python/2.7.18 mugqic/Strelka2/2.9.10 && \
rm -r -f pairedVariants/strelka2_germline && \
python $STRELKA2_HOME/bin/configureStrelkaGermlineWorkflow.py \
  --bam alignment/HCC1395BL_normal/HCC1395BL_normal.subset.bam \
  --bam alignment/HCC1395_tumor/HCC1395_tumor.subset.bam \
  --referenceFasta ${REF}/genome/Homo_sapiens.GRCh38.fa \
  --callRegions regions.bed.gz \
  --runDir pairedVariants/strelka2_germline && \
python pairedVariants/strelka2_germline/runWorkflow.py \
  -m local  \
  -j 3 \
  -g 6 \
  --quiet

# ===== Block 20/32 =====
# Filter germline
module purge && \
module load mugqic/htslib/1.14 mugqic/vt/0.57 mugqic/bcftools/1.15 && \
zcat pairedVariants/strelka2_germline/results/variants/variants.vcf.gz  | \
sed 's/TUMOR/HCC1395_NS_T_1/g'   | \
sed 's/NORMAL/HCC1395BL_NS_N_1/g'   | \
bgzip -cf  \
  > pairedVariants/HCC1395.strelka2.germline.vcf.gz && \
tabix -pvcf pairedVariants/HCC1395.strelka2.germline.vcf.gz && \
bcftools \
  view -f PASS -Oz \
 -o pairedVariants/HCC1395.strelka2.germline.loh.vcf.gz \
 pairedVariants/HCC1395.strelka2.germline.vcf.gz

# ===== Block 21/32 =====
# Decompose and normalize 4 somatic vcf files (AD field is Number=R in the VCF spec but use Number=., so we need to fix that first)
for i in pairedVariants/HCC1395.*.somatic.vcf.gz; do \
OUT=$(echo "$i" | sed 's#somatic#somatic.vt#g') ; \
echo $i ; \
module purge && \
module load mugqic/htslib/1.14 mugqic/vt/0.57 && \
zcat $i | \
sed 's/ID=AD,Number=./ID=AD,Number=R/' | \
vt decompose -s - | \
vt normalize \
    -r ${REF}/genome/Homo_sapiens.GRCh38.fa \
    - | \
bgzip -cf > ${OUT} && \
tabix -p vcf ${OUT} ; \
done

# ===== Block 22/32 =====
# Decompose and normalize 3 germline vcf files 
for i in pairedVariants/HCC1395.*.germline.loh.vcf.gz; do \
OUT=$(echo "$i" | sed 's#germline.loh#germline.loh.vt#g') ; \
echo $i ; \
module purge && \
module load mugqic/htslib/1.14 mugqic/vt/0.57 && \
zcat $i | \
sed 's/ID=AD,Number=./ID=AD,Number=R/' | \
vt decompose -s - | \
vt normalize \
    -r ${REF}/genome/Homo_sapiens.GRCh38.fa \
    - | \
bgzip -cf > ${OUT} && \
tabix -p vcf ${OUT} ; \
done

# ===== Block 23/32 =====
# Unified callset (order matters, INFO and FORMAT fields are added in the order of the input files, keep variants found by 2 callers)
module purge && \
module load mugqic/bcbio.variation.recall/0.2.6 mugqic/bcftools/1.15 mugqic/java/openjdk-jdk1.8.0_72 && \
$BCBIO_VARIATION_RECALL_HOME/bcbio.variation.recall ensemble \
  --cores 2 --numpass 1 \
  --names mutect2,strelka2,vardict,varscan2 \
  pairedVariants/HCC1395.ensemble.somatic.vcf.gz \
  ${REF}/genome/Homo_sapiens.GRCh38.fa \
  pairedVariants/HCC1395.mutect2.somatic.vt.vcf.gz \
  pairedVariants/HCC1395.strelka2.somatic.vt.vcf.gz \
  pairedVariants/HCC1395.vardict.somatic.vt.vcf.gz \
  pairedVariants/HCC1395.varscan2.somatic.vt.vcf.gz && \
zgrep -v "^#" pairedVariants/HCC1395.ensemble.somatic.vcf.gz | wc -l

# ===== Block 24/32 =====
# Unified callset (order matters, INFO and FORMAT fields are added in the order of the input files, keep variants found by 2 callers)
module purge && \
module load mugqic/bcbio.variation.recall/0.2.6 mugqic/bcftools/1.15 mugqic/java/openjdk-jdk1.8.0_72 && \
$BCBIO_VARIATION_RECALL_HOME/bcbio.variation.recall ensemble \
  --cores 2 --numpass 1 \
  --names strelka2,vardict,varscan2 \
  pairedVariants/HCC1395.ensemble.germline.vcf.gz \
  ${REF}/genome/Homo_sapiens.GRCh38.fa \
  pairedVariants/HCC1395.strelka2.germline.loh.vt.vcf.gz \
  pairedVariants/HCC1395.vardict.germline.loh.vt.vcf.gz \
  pairedVariants/HCC1395.varscan2.germline.loh.vt.vcf.gz && \
zgrep -v "^#" pairedVariants/HCC1395.ensemble.germline.vcf.gz | wc -l

# ===== Block 25/32 =====
# CPSR vcf prep (retain calls from all three caller with tumor and normal depth >=10 and tumor and normal variant allele frequency >=0.05)
module purge && \
module load mugqic/mugqic_tools/2.12.7 mugqic/python/3.10.4 mugqic/bcftools/1.15 mugqic/htslib/1.14 && \
python3 $PYTHON_TOOLS/format2pcgr.py \
  -i pairedVariants/HCC1395.ensemble.germline.vcf.gz \
  -o pairedVariants/HCC1395.ensemble.germline.format.vcf.gz \
  -f 3 \
  -v germline \
  -t HCC1395_NS_T_1 && \
bcftools \
  view -Oz -i'TDP>=10 && TVAF>=0.05 && NDP>=10 && NVAF>=0.05' \
  pairedVariants/HCC1395.ensemble.germline.format.vcf.gz | \
bcftools \
  view -Oz -s ^HCC1395_NS_T_1 | \
bcftools \
  sort -Oz \
 -o pairedVariants/HCC1395.ensemble.germline.flt.vcf.gz && \
tabix -pvcf  \
    pairedVariants/HCC1395.ensemble.germline.flt.vcf.gz

# ===== Block 26/32 =====
module purge && \
module load mugqic/pcgr/2.3.1 && \
mkdir -p pairedVariants/cpsr && \
cpsr --force_overwrite --secondary_findings --gwas_findings --pgx_findings --panel_id 0 \
    --input_vcf pairedVariants/HCC1395.ensemble.germline.flt.vcf.gz \
    --refdata_dir $PCGR_DATA \
    --vep_dir $PCGR_VEP_CACHE \
    --output_dir pairedVariants/cpsr \
    --genome_assembly grch38 \
    --sample_id HCC1395

# ===== Block 27/32 =====
# PCGR vcf prep (retain calls from all two caller with tumor and normal depth >=10 and tumor and normal variant allele frequency >=0.05)
module purge && \
module load mugqic/mugqic_tools/2.12.7 mugqic/python/3.10.4 mugqic/bcftools/1.15 mugqic/htslib/1.14 && \
python3 $PYTHON_TOOLS/format2pcgr.py \
        -i pairedVariants/HCC1395.ensemble.somatic.vcf.gz \
        -o pairedVariants/HCC1395.ensemble.somatic.format.vcf.gz \
        -f 2 \
        -v somatic \
        -t HCC1395_NS_T_1  && \
bcftools \
  view -Oz -i'TDP>=10 && TVAF>=0.05 && NDP>=10 && NVAF<=0.05' \
  -o pairedVariants/HCC1395.ensemble.somatic.flt.vcf.gz \
 pairedVariants/HCC1395.ensemble.somatic.format.vcf.gz && \
tabix -pvcf  \
  pairedVariants/HCC1395.ensemble.somatic.flt.vcf.gz 

# ===== Block 28/32 =====
## PCGR command (-tumor_site 6 is for breast cancer)
module purge && \
module load mugqic/pcgr/2.3.1 && \
mkdir -p pairedVariants/pcgr && \
pcgr --force_overwrite \
    --tumor_site 6 \
    --assay WGS \
    --call_conf_tag TAL --tumor_dp_tag TDP --tumor_af_tag TVAF --tumor_dp_min 10 --tumor_af_min 0.05 \
    --control_dp_tag NDP --control_af_tag NVAF --control_dp_min 10 --control_af_max 0.05 \
    --input_vcf pairedVariants/HCC1395.ensemble.somatic.flt.vcf.gz  \
    --input_cpsr pairedVariants/cpsr/HCC1395.cpsr.grch38.classification.tsv.gz \
    --input_cpsr_yaml pairedVariants/cpsr/HCC1395.cpsr.grch38.conf.yaml \
    --refdata_dir $PCGR_DATA \
    --vep_dir $PCGR_VEP_CACHE \
    --output_dir pairedVariants/pcgr \
    --genome_assembly grch38 \
    --sample_id HCC1395

# ===== Block 29/32 =====
# NOTE: interactive pager, for manual inspection only. Not executed here.
# less -S pairedVariants/pcgr/HCC1395.pcgr.grch38.pass.vcf.gz

# ===== Block 30/32 =====
# bcftools (use the VEP plugin to extract the CSQ field)
 module purge && \
 module load mugqic/bcftools/1.23 && \
 bcftools +split-vep -r "13:32339132" \
 pairedVariants/pcgr/HCC1395.pcgr.grch38.pass.vcf.gz \
 -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Feature\t%Consequence\t%HGVSc\t%HGVSp\t%CANONICAL\t%MANE_SELECT\n' \
 -d \
 -A tab \
 2>/dev/null

# ===== Block 31/32 =====
# Coverage Track
module purge
module load mugqic/igvtools/2.3.14

for bam in alignment/HCC1395*/*subset.bam
do
    igvtools count \
        -f min,max,mean \
        "${bam}" \
        "${bam}.tdf" \
        ${REF}/genome/Homo_sapiens.GRCh38.fa
done

# ===== Block 32/32 =====
exit
