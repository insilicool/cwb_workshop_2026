Because a single variant can carry a dozen different transcript-level consequences, PCGR (via VEP) needs to choose one to lead with in the report. 
You can see this decision recorded directly in the INFO field for our transcript of interest:

`CANONICAL=YES` — this is flagged as the Ensembl canonical transcript for BRCA2
`MANE_SELECT=NM_000059.4` — it's also the MANE Select transcript, the single RefSeq/Ensembl-agreed reference transcript recommended for clinical reporting across the genome
`PICK=1` — VEP's own --pick algorithm has explicitly selected this transcript's consequence as the one and only representative annotation, using a ranked set of criteria (canonical status, biotype, APPRIS annotation, transcript support level, and so on)

PCGR then takes this "picked" transcript's values and flattens them out of the crowded CSQ string into their own individual, easy-to-parse INFO tags 
— you can see this later in the same INFO field: 

```
Consequence=stop_gained;IMPACT=HIGH;SYMBOL=BRCA2;...;HGVSc=ENST00000380152.8:c.4777G>T;HGVSp=ENSP00000369497.3:p.Glu1593Ter;CANONICAL=YES;MANE_SELECT=NM_000059.4. 
```

This is a convenience PCGR adds on top of raw VEP output specifically so downstream reporting code (and you, reading the file) doesn't have to parse the full multi-transcript CSQ blob just to get "the" answer.
