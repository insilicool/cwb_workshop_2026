There are 14 transcripts associated with this BCRA2 variant. 

```
13	32339132	G	T	BRCA2	ENST00000380152	stop_gained	ENST00000380152.8:c.4777G>T	ENSP00000369497.3:p.Glu1593Ter	YES	NM_000059.4
13	32339132	G	T	BRCA2	ENST00000470094	stop_gained&NMD_transcript_variant	ENST00000470094.2:c.4777G>T	ENSP00000434898.2:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000528762	stop_gained&NMD_transcript_variant	ENST00000528762.2:c.4777G>T	ENSP00000433168.2:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000530893	stop_gained	ENST00000530893.7:c.4408G>T	ENSP00000499438.2:p.Glu1470Ter	.	.
13	32339132	G	T	BRCA2	ENST00000544455	stop_gained	ENST00000544455.6:c.4777G>T	ENSP00000439902.1:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000614259	stop_gained&NMD_transcript_variant	ENST00000614259.2:c.4777G>T	ENSP00000506251.1:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000665585	stop_gained&NMD_transcript_variant	ENST00000665585.2:c.4777G>T	ENSP00000499570.2:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000666593	stop_gained&NMD_transcript_variant	ENST00000666593.2:c.4777G>T	ENSP00000499256.2:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000680887	stop_gained	ENST00000680887.1:c.4777G>T	ENSP00000505508.1:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000700202	stop_gained	ENST00000700202.2:c.4777G>T	ENSP00000514856.2:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000713677	3_prime_UTR_variant&NMD_transcript_variant	ENST00000713677.1:c.*4416G>T	.	.	.
13	32339132	G	T	BRCA2	ENST00000713678	stop_gained	ENST00000713678.1:c.4777G>T	ENSP00000518981.1:p.Glu1593Ter	.	.
13	32339132	G	T	BRCA2	ENST00000713679	stop_gained&NMD_transcript_variant	ENST00000713679.1:c.4777G>T	ENSP00000518982.1:p.Glu1593Ter	.	.

```

The key point here is: the same variant, different consequences on different transcripts

This is where VEP's per-transcript annotation really matters, and this variant is a great example of it.

- ENST00000470094 (biotype = nonsense_mediated_decay) → consequence is stop_gained&NMD_transcript_variant — the transcript itself is flagged as a nonsense-mediated-decay target, 
so this "stop gain" sits on a transcript that's biologically expected to be degraded, not translated
- ENST00000713677 (also NMD biotype) → consequence is 3_prime_UTR_variant&NMD_transcript_variant — on this transcript, the same genomic position doesn't even fall in the coding sequence, 
because the exon/intron structure is different, and gets called something far less severe

So the exact same genomic change (13:32339132 G>T) is simultaneously a HIGH-impact stop-gain on some transcripts and a MODIFIER-impact 3′ UTR variant on another, 
purely because of how each transcript is structured. This is exactly why transcript choice matters when reporting a variant — reporting "3′ UTR variant" instead of "stop_gained" 
for the same base change would completely misrepresent the biological impact.