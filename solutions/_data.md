The first thing to do is to download the data and check data integrity after download e.g md5sum

The second thing to do is to verify the quality of the data.  
Typically can be done by running FastQC on the data directly or use developed pipeline like GenPipes which uses the trimmer fastp to trim and do quality assessment.  
This will give you a good idea of the quality of the data and if there are any issues that need to be addressed before proceeding with analysis.

Next run the alignment step using a tool like BWA-mem2 to align the reads to a reference genome and run post-alignment steps: minimally mark duplicate

Again assess the quality of the alignment using tools like mosdepth or samtools flagstat to check for mapping quality, coverage, and other metrics.

Other quality metrics more specfic to cancer which can impact variant assessment include:
Concordance of normal and tumor samples, and contamination, using tools like [Conpair](https://github.com/nygenome/Conpair)
Ploidy and purity estimates, using tools likes [Purples](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md) and [Sequenza](https://sequenzatools.bitbucket.io/#/home)

High contamination or low purity can have an impact on variant calling and interpretation, so it is important to assess and understand the consequences of these metrics before proceeding with downstream analysis.

Finally, after all quality checks and assessments are done, you can proceed with variant calling and other downstream analyses.