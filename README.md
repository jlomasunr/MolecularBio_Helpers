# MolecularBio_Helpers
Custom scripts to help with molecular biology tasks.

Functionality includes:
- Reference guided plasmid sequence verification
- Codon optimization for DNA synthesis

# Plasmid comparison (ComPlasmid.py)
ComPlasmid helps one verify a plasmid sequence (i.e. from a whole plasmid sequencing service) against a target plasmid (i.e. from *in silico* cloning). The task is complicated by the fact that a sequenced and assembled plasmid may have a different starting location (index) and orientation than the reference, and may contain sequencing errors. One may also be more interested in specific features, therefore annotations from the reference need to be applied to the query.

**Re-indexing:**
A stringent local alignment of the query and reference plasmid is used to re-index the query to match orientation and start location of the reference as closely as possible. Both forward and reverse orientations of the query are searched and the orientation with the highest scoring local alignment is chosen.

**Feature transfer**: 
If the reference is provided in GenBank format, features are extracted and locally aligned to the re-indexed query. The bounds of the local alignments on the query are taken as the new coordinates of the annotation.
## Environemnt
```
conda install -c bioconda biopython
```

## Usage
```
ComPlasmid.py [-h] [--transfer] subject_file query_file

positional arguments:
  subject_file  Subject file (fasta or genbank)
  query_file    Query file (fasta or genbank)

options:
  -h, --help    show this help message and exit
  --transfer    Map subject features to query file
```

***Input:***
Fasta or GenBank reference and query plasmids. Features can only be mapped if a GenBank reference file is provided.

***Output:***
Re-indexed query plasmid in fasta and genbank format. A report of the re-indexing and feature mapping results.

### Example
```
ComPlasmid.py --transfer data/AvBa5-Ats_subject.gbk data/AvBa5Ats_query.gbk
```

----
----

# Codon Usage Tables
Arabidopsis codon usage proportions (https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3702&aa=1&style=N). Values represent fractions among synonymous codons

Glycine max codon usage proportions (https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3847&aa=1&style=N). Values represent fractions among synonymous codons