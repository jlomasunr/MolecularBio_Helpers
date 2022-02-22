#!/usr/bin/env python

# Usage: ./codonOptimize.py [fasta_file]
# Script to codon optimize a fasta file for Arabidopsis thaliana. Codons in the
# input file are replaced with the most frequently used Arabidopsis codons.

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Arabidopsis codon usage proportions (https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=3702&aa=1&style=N)
# Values represent fractions among synonymous codons
AT_CodonUsage = {
    "TTT": 0.51,"TTC": 0.49,"TTA": 0.14,"TTG": 0.22,
    "CTT": 0.26,"CTC": 0.17,"CTA": 0.11,"CTG": 0.11,
    "ATT": 0.41,"ATC": 0.35,"ATA": 0.24,"ATG": 1.00,
    "GTT": 0.40,"GTC": 0.19,"GTA": 0.15,"GTG": 0.26,
    "TCT": 0.28,"TCC": 0.13,"TCA": 0.20,"TCG": 0.10,
    "AGT": 0.16,"AGC": 0.13,"CCT": 0.38,"CCC": 0.11,
    "CCA": 0.33,"CCG": 0.18,"ACT": 0.34,"ACC": 0.20,
    "ACA": 0.31,"ACG": 0.15,"GCT": 0.43,"GCC": 0.16,
    "GCA": 0.27,"GCG": 0.14,"TAT": 0.52,"TAC": 0.48,
    "CAT": 0.61,"CAC": 0.39,"CAA": 0.56,"CAG": 0.44,
    "AAT": 0.52,"AAC": 0.48,"AAA": 0.49,"AAG": 0.51,
    "GAT": 0.68,"GAC": 0.32,"GAA": 0.52,"GAG": 0.48,
    "GAA": 0.52,"GAG": 0.48,"TGT": 0.60,"TGC": 0.40,
    "TGG": 1.00,"CGT": 0.17,"CGC": 0.07,"CGA": 0.12,
    "CGG": 0.09,"AGA": 0.35,"AGG": 0.20,"GGT": 0.34,
    "GGC": 0.14,"GGA": 0.37,"GGG": 0.16,"TAA": 0.36,
    "TAG": 0.20,"TGA": 0.44}
CodonsDict = {
    "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
    "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
    "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
    "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
    "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
    "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
    "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
    "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
    "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
    "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
    "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
    "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
    "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
    "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
    "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
    "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0
    }

SynonymousCodons = {
    "CYS": ["TGT", "TGC"],
    "ASP": ["GAT", "GAC"],
    "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "GLN": ["CAA", "CAG"],
    "MET": ["ATG"],
    "ASN": ["AAC", "AAT"],
    "PRO": ["CCT", "CCG", "CCA", "CCC"],
    "LYS": ["AAG", "AAA"],
    "STOP": ["TAG", "TGA", "TAA"],
    "THR": ["ACC", "ACA", "ACG", "ACT"],
    "PHE": ["TTT", "TTC"],
    "ALA": ["GCA", "GCC", "GCG", "GCT"],
    "GLY": ["GGT", "GGG", "GGA", "GGC"],
    "ILE": ["ATC", "ATA", "ATT"],
    "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "HIS": ["CAT", "CAC"],
    "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "TRP": ["TGG"],
    "VAL": ["GTA", "GTC", "GTG", "GTT"],
    "GLU": ["GAG", "GAA"],
    "TYR": ["TAT", "TAC"]
}

def findAmino(codon):
    for amino in SynonymousCodons:
        if codon in SynonymousCodons[amino]:
            return(amino)

def topCodonForAmino(amino):
    topAT = [SynonymousCodons[amino][0], AT_CodonUsage[SynonymousCodons[amino][0]]]
    for codon in SynonymousCodons[amino]:
        if topAT[1] < AT_CodonUsage[codon]:
            topAT =  [codon, AT_CodonUsage[codon]]
    return(topAT[0])

def optimizeCodons(record):
    sequence = str(record.seq)
    desc = "Codon optimized for Arabidopsis thaliana"
    optimized_sequence = ""

    for i in range(0, len(sequence), 3):
        codon = sequence[i : i + 3]
        if codon in CodonsDict.keys():
            amino = findAmino(codon)
            optimized_sequence = optimized_sequence + topCodonForAmino(amino)
        else:
            print("Could not find codon!....Exiting")
            return()

    optimized_record = SeqRecord(Seq(optimized_sequence),
                                 id = record.id,
                                 description = desc)

    return(optimized_record)

def main():
    FILE = sys.argv[1]
    fasta = SeqIO.to_dict(SeqIO.parse(FILE, "fasta"))
    optimized_records = []

    for i in range(0,len(fasta.keys())):
        record = fasta[list(fasta)[i]]
        optimized = optimizeCodons(record)
        print(optimized)
        optimized_records.append(optimized)

    SeqIO.write(optimized_records, "optimized.fa", "fasta")

main()
