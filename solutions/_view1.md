Using samtools view -H we can see the header of the bam file. 

The header contains information about the reference genome, read groups, and other metadata.

@HD VN:1.6  SO:coordinate tells us that the file is in SAM format version 1.6 and that the alignments are sorted by coordinate.

@SQ SN:chr1 LN:248956422 tells us that the reference genome (GRCh38) contains a sequence named chr1 with a length of 248,956,422 base pairs.

@PG Describes what programs, version and command line were used to generate the BAM file.
- In this case BWA-mem2 v2.2.1was used to align the trimmed (fastp) reads to the reference genome.
- GATK MarkDuplicates v4.6.0.0 was used to mark duplicate reads in the BAM file.
- Samtools view using regions.snv.bed was used to subset the BAM file to only include reads that overlap the regions specified in the BED file.
  Caveat: Subsetting removes reads whose mate maps outside the region, which can subtly affect callers that use mate-pair/insert-size information (like Strelka2 and Mutect2 both do some local reassembly)

A full description of the flags can be found in the SAM specification
http://samtools.sourceforge.net/SAM1.pdf
