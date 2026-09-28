Description: Run VarScan2 to call somatic variants from mpileup file

pairedVariants/HCC1395.varscan2.mpileup  #input mpileup file
pairedVariants/HCC1395.varscan2          #output prefix
--min-coverage 3                         #minmium read depth (ref + alt; default 8)
--min-var-freq 0.05                      #minimum variant allele frequency (default 0.20)
--p-value 0.10                           #p-value threshold for calling variants (default 0.99)
--somatic-p-value 0.05                   #p-value threshold for calling somatic variants (default 0.05)
--strand-filter 0                        #disable filtering of variants with >90% support on one strand (default 1)
--output-vcf 1                           #output VCF format (default 0)
--mpileup 1                              #input file is in mpileup format (default 0)

bcftools concat -a                       #merge snp and indel VCFs into a single VCF  

sed 's/TUMOR/HCC1395_NS_T_1/g'           #Rename generic TUMOR and NORMAL sample names in VCF to something more meaningful
sed 's/NORMAL/HCC1395BL_NS_N_1/g' 