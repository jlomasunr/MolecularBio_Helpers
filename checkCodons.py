#!/usr/bin/env python

# Usage: ./checkCodons.py -s [species] [original_fasta] [optimized_fasta]
desc = '''
Script to check a codon optimized fasta file against the original sequences.
'''

import sys
import pandas as pd
import dnachisel
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import argparse
import os
import json

arg_parser = argparse.ArgumentParser(description=desc)
arg_parser.add_argument('original',
                        type=str,
                        help='Path to original fasta file')
arg_parser.add_argument('optimized',
                        type=str,
                        help='Path to optimized fasta file')
arg_parser.add_argument('-s', 
                        metavar="species", 
                        type=str, 
                        help="Species codon usage to use ['A_thaliana', 'G_max']", 
                        default="A_thaliana", 
                        choices=['A_thaliana', 'G_max'])

args = arg_parser.parse_args()

if args.s == "G_max":
    codonUsageFile = "G_max_codonUsage.json"
else:
    codonUsageFile = "A_thaliana_codonUsage.json"

scriptDir = os.path.realpath(__file__).replace("checkCodons.py", "")
with open(scriptDir + codonUsageFile) as cUF:
    codonUsageTable = json.load(cUF)

codonUsage = {}
for aa in codonUsageTable.keys():
    for codon in codonUsageTable[aa].keys():
        codonUsage[codon] = codonUsageTable[aa][codon]

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

def countCodons(dna_sequence):
    tempCodonCounts = CodonsDict
    tempCodonFreq = CodonsDict
    # Compute codon counts
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i : i + 3]
        if codon in tempCodonCounts.keys():
            tempCodonCounts[codon] += 1

    # Compute frequencies
    for codon in tempCodonFreq:
        amino = findAmino(codon)
        synonymousCounts = [tempCodonCounts[c] for c in SynonymousCodons[amino]]
        if sum(synonymousCounts) == 0:
            tempCodonFreq[codon] = 0
        else:
            tempCodonFreq[codon] = tempCodonCounts[codon]/sum(synonymousCounts)

    return(tempCodonFreq)

def rankCodon(amino, codon):
    ranking = {}
    for AT_codon in SynonymousCodons[amino]:
        ranking[AT_codon] = codonUsage[AT_codon]
    ranking = pd.DataFrame.from_dict(ranking, orient='index').sort_values(0, ascending=False).reset_index()
    return(ranking.index[ranking["index"]==codon].to_list()[0] + 1)

def main():
    FILE_ORIGINAL = args.original
    FILE_OPTIMIZED = args.optimized

    original_records = SeqIO.to_dict(SeqIO.parse(FILE_ORIGINAL, "fasta"))
    optimized_records = SeqIO.to_dict(SeqIO.parse(FILE_OPTIMIZED, "fasta"))

    #rt = {"AT_freq":codonUsage}
    rt = pd.DataFrame.from_dict([codonUsage]).transpose().sort_index()
    rt.columns = ["AT"]


    for i in range(0,len(original_records.keys())):
        original_dna = original_records[list(original_records)[i]]
        optimized_dna = optimized_records[list(optimized_records)[i]]
        try:
            dnaDiffs = dnachisel.biotools.annotate_differences(optimized_dna, original_dna)
            print(dnaDiffs)
        except:
            print("No differences")

        #alignments = pairwise2.align.globalxx(original_dna.seq, optimized_dna.seq)
        #print(format_alignment(*alignments[0]))
        counts = countCodons(str(optimized_dna.seq))

        for amino in SynonymousCodons:
            topAT = [SynonymousCodons[amino][0], codonUsage[SynonymousCodons[amino][0]]]
            topOther = [SynonymousCodons[amino][0], counts[SynonymousCodons[amino][0]], rankCodon(amino, SynonymousCodons[amino][0]), len(SynonymousCodons[amino])]
            for codon in SynonymousCodons[amino]:
                if topAT[1] < codonUsage[codon]:
                    topAT =  [codon, codonUsage[codon]]
                if topOther[1] < counts[codon]:
                    topOther = [codon, counts[codon], rankCodon(amino, codon), len(SynonymousCodons[amino])]
            if topAT[0] == topOther[0]:
                flag = ""
            else:
                flag = "---Different top codons---  "
            print("%s%s:   AT:%s   CDS:%s" % (flag, amino, topAT, topOther))

        counts = pd.DataFrame.from_dict([counts]).transpose().sort_index()
        rt[original_dna.id] = counts[0]

        original_protein = SeqRecord(original_dna.seq.translate())
        optimized_protein = SeqRecord(optimized_dna.seq.translate())
        #proteinDiffs = dnachisel.biotools.annotate_differences(optimized_protein, original_protein)
        #print(proteinDiffs)
        alignments = pairwise2.align.globalxx(original_protein.seq, optimized_protein.seq)
        print(format_alignment(*alignments[0]))

    rt = rt.reset_index()
    rt = rt.rename(columns = {'index':'codon'})
    return(rt)

main()
