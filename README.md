# MolecularBio_Helpers
Custom scripts to help with molecular biology tasks.

Functionality includes:
- [Reference guided plasmid sequence verification](#plasmid-comparison-(ComPlasmid.py))
- [Codon optimization for DNA synthesis](#codon-optimization-for-dna-synthesis)

# Plasmid comparison (ComPlasmid.py)
ComPlasmid helps one verify a plasmid sequence (i.e. from a whole plasmid sequencing service) against a target plasmid (i.e. from *in silico* cloning). The task is complicated by the fact that a sequenced and assembled plasmid may have a different starting location (index) and orientation than the reference, and may contain sequencing errors. One may also be more interested in specific features, therefore annotations from the reference need to be applied to the query.

**Re-indexing:**
A stringent local alignment of the query and reference plasmid is used to re-index the query to match orientation and start location of the reference as closely as possible. Both forward and reverse orientations of the query are searched and the orientation with the highest scoring local alignment is chosen.

**Feature transfer**: 
If the reference is provided in GenBank format, features are extracted and locally aligned to the re-indexed query. The bounds of the local alignments on the query are taken as the new coordinates of the annotation.
## Environemnt
```
mamba install -c bioconda -c conda-forge parasail-python dnachisel biopython
```

## Usage
```
ComPlasmid.py [-h] [options] subject_file query_file

positional arguments:
  subject_file  Subject file (fasta or genbank)
  query_file    Query file (fasta or genbank)

options:
  -h, --help    show this help message and exit
  --transfer    Map subject features to query file
  --strand      [forward|reverse] Only search one strand for local alignments
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
# Codon optimization for DNA synthesis
A set of scripts reliant on [BioPython](https://biopython.org/) and [DNAChisel](https://edinburgh-genome-foundry.github.io/DnaChisel/) to perfrom codon optimization with the dual objectives of efficent translation and efficent DNA synthesis. Sequences are optimized towards the fastest translating codons while avoiding common motifs that inhibit DNA synthesis.

Avoided motifs include:
- Repeated K-mers
- Homopolymers
- Hairpins

## Environment
```
conda install -c bioconda dnachisel biopython
```

## Workflow
Start with a fasta file of coding sequences to be optimized and run ```codonOptimizeSynthesis.py```. Adjust the weights and thresholds as necessary to acheive the desired result.

```
usage: codonOptimizeSynthesis.py [options] input_file

positional arguments:
  input_file            Fasta file to be optimized

options:
  -h, --help            show this help message and exit
  -o outfile            Output file name
  -s species            Species codon usage to optimize towards ['A_thaliana', 'G_max']
  --rscu RSCU           JSON file containing RSCU values
  --rscuGenerate CDS.fasta
                        Generate rscu values from CDS fasta file
  --topCodonWeight TOPCODONWEIGHT
                        Adjust weight optimizing towards highest frequency codon (default: 25)
  --rareCodonThreshold RARECODONTHRESHOLD
                        Attempt to avoid codons with frequencies below this number (default: 0.1)
  --rareCodonWeight RARECODONWEIGHT
                        Adjust weight for avoiding rare codons (default: 500)
  --repeatedKmerWeight REPEATEDKMERWEIGHT
                        Adjust weight for avoiding repeated kmers (default: 100)
  --repeatedHomopolymerWeight REPEATEDHOMOPOLYMERWEIGHT
                        Adjust weight for avoiding repeated homopolymers (default: 50)
  --hairpinWeight HAIRPINWEIGHT
                        Adjust weight for avoiding hairpins (default: 1000)
  --uniquifyKmersWeight UNIQUIFYKMERSWEIGHT
                        Adjust weight for avoiding kmer repeats of 8 bases and longer (default: 1000)
```

Next, you can verify your optimized sequences using ```checkCodons.py```. This will print out a report showing the original and optimized protein alignment and the usage of top codons over the coding sequence.

```
checkCodons.py [-h] [-s species] original optimized


positional arguments:
  original    Original fasta file
  optimized   Optimized fasta file

options:
  -h, --help  show this help message and exit
  -s species  Species codon usage to use ['A_thaliana', 'G_max']
```

## Codon Usage Tables

Codon usage tables for *A. thaliana* and *Glycine max* (soybean) are provided in the ```data/``` folder. The values represent fractions among synonymous codons. You can provide your own codon usage table via the ```--rscu``` argument.

**Note:** Many measures of codon translation efficiency exist and there may be a more effective formulation for your use case... These include (CAI, ENC, tAI, expression normalized tAI, etc.)

**Arabidopsis:** 
https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3702&aa=1&style=N

**Soybean** https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3847&aa=1&style=N