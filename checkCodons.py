#!/usr/bin/env python

# Usage: ./checkCodons.py [original_fasta] [optimized_fasta]

import sys
import pandas as pd
import dnachisel
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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

def main():
    FILE_ORIGINAL = sys.argv[1]
    FILE_OPTIMIZED = sys.argv[2]

    original_records = SeqIO.to_dict(SeqIO.parse(FILE_ORIGINAL, "fasta"))
    optimized_records = SeqIO.to_dict(SeqIO.parse(FILE_OPTIMIZED, "fasta"))

    #rt = {"AT_freq":AT_CodonUsage}
    rt = pd.DataFrame.from_dict([AT_CodonUsage]).transpose().sort_index()
    rt.columns = ["AT"]


    for i in range(0,len(original_records.keys())):
        original_dna = original_records[list(original_records)[i]]
        optimized_dna = optimized_records[list(optimized_records)[i]]
        dnaDiffs = dnachisel.biotools.annotate_differences(optimized_dna, original_dna)
        print(dnaDiffs)

        #alignments = pairwise2.align.globalxx(original_dna.seq, optimized_dna.seq)
        #print(format_alignment(*alignments[0]))
        counts = countCodons(str(optimized_dna.seq))

        for amino in SynonymousCodons:
            topAT = [SynonymousCodons[amino][0], AT_CodonUsage[SynonymousCodons[amino][0]]]
            topOther = [SynonymousCodons[amino][0], counts[SynonymousCodons[amino][0]]]
            for codon in SynonymousCodons[amino]:
                if topAT[1] < AT_CodonUsage[codon]:
                    topAT =  [codon, AT_CodonUsage[codon]]
                if topOther[1] < counts[codon]:
                    topOther = [codon, counts[codon]]
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
