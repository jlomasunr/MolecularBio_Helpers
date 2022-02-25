#!/usr/bin/env python

# Usage: ./codonOptimizeSynthesis.py [fasta_file]
# Script to codon optimize a fasta file for Arabidopsis thaliana. Codons in the
# input file are replaced with the most frequently used Arabidopsis codons.

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dnachisel import *

# Arabidopsis codon usage proportions (https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3702&aa=1&style=N)
# Values represent fractions among synonymous codons
AT_CodonUsage = {
    "F":{"TTT":0.51, "TTC": 0.49},
    "L":{"TTA": 0.14, "TTG": 0.22, "CTT": 0.26, "CTC": 0.17, "CTA": 0.11, "CTG": 0.11},
    "S":{"TCT": 0.28, "TCC": 0.13, "TCA": 0.2, "TCG": 0.1, "AGT": 0.16, "AGC": 0.13},
    "P":{"CCT": 0.38, "CCC": 0.11, "CCA": 0.33, "CCG": 0.18},
    "R":{"CGT": 0.17, "CGC": 0.07, "CGA": 0.12, "CGG": 0.09, "AGA": 0.35, "AGG": 0.2},
    "V":{"GTA": 0.15, "GTG": 0.26, "GTT": 0.4, "GTC": 0.19},
    "T":{"ACT": 0.34, "ACC": 0.2, "ACA": 0.31, "ACG": 0.15},
    "A":{"GCT": 0.43, "GCC": 0.16, "GCA": 0.27, "GCG": 0.14},
    "I":{"ATT": 0.41, "ATC": 0.35, "ATA": 0.24},
    "M":{"ATG": 1},
    "G":{"GGG": 0.16, "GGT": 0.34, "GGC": 0.14, "GGA": 0.37},
    "W":{"TGG": 1},
    "C":{"TGT": 0.6, "TGC": 0.4},
    "Y":{"TAT": 0.52, "TAC": 0.48},
    "H":{"CAT": 0.61, "CAC": 0.39},
    "Q":{"CAA": 0.56, "CAG": 0.44},
    "N":{"AAT": 0.52, "AAC": 0.48},
    "K":{"AAA": 0.49, "AAG": 0.51},
    "D":{"GAT": 0.68, "GAC": 0.32},
    "E":{"GAA": 0.52, "GAG": 0.48},
    "*":{"TAA": 0.36, "TAG": 0.2, "TGA": 0.44}
}

FILE = sys.argv[1]
fasta = SeqIO.to_dict(SeqIO.parse(FILE, "fasta"))
optimized_records = []

for i in range(0,len(fasta.keys())):
    record = fasta[list(fasta)[i]]
    print("--------------  %s --------------" % (record.id))
    # DEFINE THE OPTIMIZATION PROBLEM
    problem = DnaOptimizationProblem(
        sequence=str(record.seq),
        constraints=[
            AvoidPattern(RepeatedKmerPattern(4, 3)),
            AvoidHairpins(stem_size=9),
            EnforceTranslation()
            ],
        objectives=[CodonOptimize(original_codon_usage_table=AT_CodonUsage, codon_usage_table=AT_CodonUsage, method="harmonize_rca")]
    )
    try:
        # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE
        problem.max_random_iters = 20000
        problem.resolve_constraints()
        problem.optimize()

        # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS
        print(problem.constraints_text_summary())
        print(problem.objectives_text_summary())

        # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)
        final_sequence = problem.sequence  # string
        final_record = problem.to_record(with_sequence_edits=True)
        final_record.id = record.id

        optimized_records.append(final_record)
    except:
        print("No solution found!")

SeqIO.write(optimized_records, "optimizedSynthesis.fa", "fasta")




