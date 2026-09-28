*For vardictjava*
Filter thresholds for variant calling:
-f 0.03 : allele frequency threshold default: 0.05
-Q 10   : minimum mapping quality default: 20
-q 25   : minimum base quality default: 25

Layout for columns in output file to pipe into testsomatic.R and var2vcf_paired.pl
-c 1  : chromosome column 1 
-S 2  : Start position column 2
-E 3  : End position column 3
-g 4  : The column for gene name, or segment annotation 
-th 3 : number of threads to use default: 1

*For var2vcf_paired.pl*
-f 0.03 : allele frequency threshold default: 0.02
-P 0.9  : The maximum p-value.  Default to 0.05
-m 4.25 : The maximum mean mismatches allowed.  Default: 5.25, or if a variant is supported by reads with more than 5.25 mismathes, it'll be considered false positive.  Mismatches don't includes indels in the alignment.
-M      : If set, output only candidate somatic
