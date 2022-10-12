#!/usr/bin/env python

desc='''
Script to codon optimize a fasta file for Arabidopsis thaliana. The program optimizes 
DNA sequences to balance the reduction of sequence complexity for DNA synthesis and 
top codon usage for improved translation. Input: DNA sequences in fasta format
'''

# Improvements:
# Print summary of non-passing constraints and their location
# Compute sequence complexity score...?

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dnachisel import *
import argparse
import json
import logging 
import os
from rich.logging import Console, RichHandler
from codonUsageCalculations import calculateCAI
from codonUsageCalculations import calculateRSCU

logging.basicConfig(level=logging.INFO,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=Console(stderr=True))])

arg_parser = argparse.ArgumentParser(description=desc)
arg_parser.add_argument('file',
                        type=str,
                        help='Path to fasta file to be optimized')
arg_parser.add_argument('-o',
                        metavar='outfile',
                        type=str,
                        default='optimized.fa',
                        help='Output file name')
arg_parser.add_argument('-s', 
                        metavar="species", 
                        type=str, 
                        help="Species codon usage to use ['A_thaliana', 'G_max']", 
                        default="A_thaliana", 
                        choices=['A_thaliana', 'G_max'])
arg_parser.add_argument('--topCodonWeight',
                        type=int,
                        default=25,
                        help='Adjust weight optimizing towards highest frequency codon')
arg_parser.add_argument('--rareCodonThreshold',
                        type=float,
                        default=0.1,
                        help='Attempt to avoid codons with frequencies below this number')
arg_parser.add_argument('--rareCodonWeight',
                        type=int,
                        default=500,
                        help='Adjust weight for avoiding rare codons')
arg_parser.add_argument('--repeatedKmerWeight',
                        type=int,
                        default=100,
                        help='Adjust weight for avoiding repeated kmers')
arg_parser.add_argument('--repeatedHomopolymerWeight',
                        type=int,
                        default=50,
                        help='Adjust weight for avoiding repeated homopolymers')
arg_parser.add_argument('--hairpinWeight',
                        type=int,
                        default=1000,
                        help='Adjust weight for avoiding hairpins')
arg_parser.add_argument('--uniquifyKmersWeight',
                        type=int,
                        default=1000,
                        help='Adjust weight for avoiding kmer repeats of 8 bases and longer')
arg_parser.add_argument('--rscu',
                        type=str,
                        default=None,
                        help='JSON file containing RSCU values for species')
arg_parser.add_argument('--rscuGenerate',
                        type=str,
                        default=None,
                        help='Generate rscu values from CDS fasta file')

args = arg_parser.parse_args()

# Load the codon usage table
if args.s == "G_max":
    codonUsageFile = "G_max_codonUsage.json"
else:
    codonUsageFile = "A_thaliana_codonUsage.json"

scriptDir = os.path.realpath(__file__).replace("codonOptimizeSynthesis.py", "")

logging.info(f"Loading codon usage: {scriptDir}/{args.s}")
with open(scriptDir + codonUsageFile) as cUF:
    codonUsageTable = json.load(cUF)


# Load the fasta file
logging.info(f"Loading sequences from '{args.file}'")
fasta = SeqIO.to_dict(SeqIO.parse(args.file, "fasta"))

# Load RSCU file
if args.rscuGenerate != None:
    rscu_dict = calculateRSCU(args.rscuGenerate)
elif args.rscu != None:
    logging.info(f"Loading RSCU values from '{args.rscu}'")
    with open(args.rscu) as rscu:
        rscu_dict = json.load(rscu)

# Optimization objectives
problem_objectives = [
    CodonOptimize(codon_usage_table=codonUsageTable, method="use_best_codon", boost=args.topCodonWeight),
    AvoidRareCodons(args.rareCodonThreshold, codon_usage_table=codonUsageTable, boost=200),
    AvoidPattern(RepeatedKmerPattern(3, 3), boost=args.repeatedKmerWeight),
    AvoidPattern(RepeatedKmerPattern(2, 8), boost=args.repeatedKmerWeight),
    AvoidPattern(RepeatedKmerPattern(2, 9), boost=args.repeatedKmerWeight),
    AvoidPattern(RepeatedKmerPattern(2, 10), boost=args.repeatedKmerWeight),
    AvoidPattern(RepeatedKmerPattern(2, 11), boost=args.repeatedKmerWeight),
    AvoidPattern(RepeatedKmerPattern(2, 12), boost=args.repeatedKmerWeight),
    AvoidPattern(RepeatedKmerPattern(2, 13), boost=args.repeatedKmerWeight),
    AvoidPattern(RepeatedKmerPattern(2, 14), boost=args.repeatedKmerWeight),
    AvoidPattern(RepeatedKmerPattern(2, 15), boost=args.repeatedKmerWeight),
    AvoidPattern(HomopolymerPattern("A", 6), boost=args.repeatedHomopolymerWeight),
    AvoidPattern(HomopolymerPattern("T", 6), boost=args.repeatedHomopolymerWeight),
    AvoidPattern(HomopolymerPattern("G", 6), boost=args.repeatedHomopolymerWeight),
    AvoidPattern(HomopolymerPattern("C", 6), boost=args.repeatedHomopolymerWeight),
    AvoidPattern(HomopolymerPattern("A", 5), boost=args.repeatedHomopolymerWeight),
    AvoidPattern(HomopolymerPattern("T", 5), boost=args.repeatedHomopolymerWeight),
    AvoidPattern(HomopolymerPattern("G", 5), boost=args.repeatedHomopolymerWeight),
    AvoidPattern(HomopolymerPattern("C", 5), boost=args.repeatedHomopolymerWeight),
    AvoidHairpins(stem_size=9, boost=args.hairpinWeight),
    AvoidHairpins(stem_size=10, boost=args.hairpinWeight),
    AvoidHairpins(stem_size=11, boost=args.hairpinWeight),
    AvoidHairpins(stem_size=12, boost=args.hairpinWeight),
    UniquifyAllKmers(8, boost=args.uniquifyKmersWeight)
]

# Optimize each entry in the fasta file
optimized_records = []
for i in range(0,len(fasta.keys())):
    record = fasta[list(fasta)[i]]
    logging.info(f"Optimizing {record.id}")

    problem = DnaOptimizationProblem(
        sequence=str(record.seq),
        constraints=[EnforceTranslation()],
        objectives=problem_objectives)

    try:
        # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE
        problem.max_random_iters = 20000
        problem.optimize()
        problem.resolve_constraints(final_check=True)

        # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS
        print(problem.objectives_text_summary())
        print(problem.constraints_text_summary())

        # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)
        final_record = problem.to_record(with_sequence_edits=True)
        final_record.id = record.id
        final_record.description = ""

        print(str(final_record.seq))

        # Check RCSU values
        print(f"Original sequence CAI: {calculateCAI(record, rscu_dict)}")
        print(f"Optimized sequence CAI: {calculateCAI(final_record, rscu_dict)}")

        optimized_records.append(final_record)
    except:
        logging.error("No solution found!")

SeqIO.write(optimized_records, args.o, "fasta")

logging.info(f"Done! Optimized sequences written to {args.o}")