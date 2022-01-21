# PGPT-Pred
The PGPT-Pred tool (https://plabase.informatik.uni-tuebingen.de/pb/form.php?var=PGPT-Pred) is part of the web resource PLaBAse (https://plabase.informatik.uni-tuebingen.de/pb/plabase.php) and allows annotation of bacterial plant growth-promoting traits (proteins), short "PGPTs" of single genomes, using blastp+hmmer or IMG-KEGG-annoation Mapper against the PGPT ontology.

When applying PGPT-Pred via the PLaBAse always cite the respective reference:

Patz S, Gautam A, Becker M, Ruppel S, Rodríguez-Palenzuela P, Huson DH. PLaBAse: A comprehensive web resource for analyzing the plant growth-promoting potential of plant-associated bacteria. bioRxiv 2021, https://doi.org/10.1101/2021.12.13.472471 (preprint)


## Input formats:
1. blastp+hmmer (PGPTblhm):            genomic protein sequences in FASTA format (sorted by genomic location)
2. IMG-KEGG-annoation Mapper (PGPTpy): genomic protein KEGG annotations (received by IMG Server, or customer format)

KEGG-Costumer formats:
- Please have a look into the Manual for possible KEGG annotation formats, that are accepted: https://plabase.informatik.uni-tuebingen.de/pb/manual.php

## Pipelines applied:
### 1. Protein/KEGG Annotation/Mapping against PGPT ontology
   1. blastp+hmmer (PGPTblhm):            genomic protein sequences are aligned against proteins associated with the PGPT ontology and respective PFAM domain comparison is achieved by hmmer against the PFAM domains using pgpt_blhm.py
   2. IMG-KEGG-annoation Mapper (PGPTpy): KEGG annotations (one per protein only) are mapped against the PGPT ontology, using pgpt_comp_fun_ascii_v2.py

### 2. Pie Chart generation 
   1. based on blastp+hmmer results or all blast hits (ignoring pfam comparison) of PGPTs, by applying:
   2. based on KEGG-PGPT mapping, by applying: 

### 3. Krona Plot generation
   1. based on blastp+hmmer results or all blast hits (ignoring pfam comparison) of PGPTs, by applying:
   2. based on KEGG-PGPT mapping, by applying: 

## Results:
1. Download:   Summary file listing all blastp+hmmer or KEGG-mapped hits of PGPTs
2. Pie Chart:  Pie Chart  summarizing either all blastp+hmmer hits or all blast hits (ignoring pfam comparison) or KEGG-PGPT hits in a percentage scale on Ontology level 2
3. Krona Plot: Krona Plot giving an hierachical overview of either all blastp+hmmer hits or all blast hits (ignoring pfam comparison) of PGPTs or KEGG-PGPT hits across all hierarchical levels
