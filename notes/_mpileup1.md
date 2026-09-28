Samtools:

	-d 1000 : max per-file depth; avoids excessive memory usage [8000] ; 
	-B : disable BAQ (per-Base Alignment Quality) ; 
	-q 10 : skip alignments with mapQ smaller than 10 [0]; 
    -Q 10 : skip bases with baseQ smaller than 10 [13];


