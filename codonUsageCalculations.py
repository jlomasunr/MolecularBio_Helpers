#!/usr/bin/env python

import sys
import os
import logging
import warnings
import re
from rich.logging import Console, RichHandler
import Bio.SeqIO as SeqIO
import CAI
import json

warnings.filterwarnings("error")

logging.basicConfig(level=logging.INFO,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=Console(stderr=True))])

if len(sys.argv) < 2:
	sys.exit("Usage: calculateRSCU.py [cds.fa]")

SynonymousCodons = {
    "CYS": {"TGT":0, "TGC":0},
    "ASP": {"GAT":0, "GAC":0},
    "SER": {"TCT":0, "TCG":0, "TCA":0, "TCC":0, "AGC":0, "AGT":0},
    "GLN": {"CAA":0, "CAG":0},
    "MET": {"ATG":0},
    "ASN": {"AAC":0, "AAT":0},
    "PRO": {"CCT":0, "CCG":0, "CCA":0, "CCC":0},
    "LYS": {"AAG":0, "AAA":0},
    "STOP": {"TAG":0, "TGA":0, "TAA":0},
    "THR": {"ACC":0, "ACA":0, "ACG":0, "ACT":0},
    "PHE": {"TTT":0, "TTC":0},
    "ALA": {"GCA":0, "GCC":0, "GCG":0, "GCT":0},
    "GLY": {"GGT":0, "GGG":0, "GGA":0, "GGC":0},
    "ILE": {"ATC":0, "ATA":0, "ATT":0},
    "LEU": {"TTA":0, "TTG":0, "CTC":0, "CTT":0, "CTG":0, "CTA":0},
    "HIS": {"CAT":0, "CAC":0},
    "ARG": {"CGA":0, "CGC":0, "CGG":0, "CGT":0, "AGG":0, "AGA":0},
    "TRP": {"TGG":0},
    "VAL": {"GTA":0, "GTC":0, "GTG":0, "GTT":0},
    "GLU": {"GAG":0, "GAA":0},
    "TYR": {"TAT":0, "TAC":0}
}

def calculateRSCU(file):
	try:
		# Compute RSCU
		fasta = [str(i.seq) for i in list(SeqIO.to_dict(SeqIO.parse(file, "fasta")).values())]
		rscu = CAI.RSCU(fasta)
	except:
		logging.info("Not all sequences divisible by three...Removing incompatable sequences")
		with open("removedIDs", 'w') as to_remove:
			for record in SeqIO.parse(file, "fasta"):
				try:
					record.seq.translate(cds=True)
				except:
					logging.warning(f"Could not translate {record.id}! Removing sequence...")
					to_remove.write(f"{str(record.id)}\n")
		
		os.system(f"seqkit grep --quiet -v -n -r -f removedIDs {file} > {file}.cleaned")
		fasta = [str(i.seq) for i in list(SeqIO.to_dict(SeqIO.parse(file + ".cleaned", "fasta")).values())]
		rscu = CAI.RSCU(fasta)

	# Write RSCU to file
	cai_json = json.dumps(rscu, indent=4)
	rscu_file = sys.argv[1].replace(".fa", ".rscu")
	with open(rscu_file, "w") as r:
		r.write(cai_json)
	logging.info(f"RSCU values written to {rscu_file}")

	return(rscu)

def calculateCAI(record, rscu):
	cai = CAI.CAI(str(record.seq), RSCUs=rscu)
	return(cai)

def translateSeqToRSCU(seq, rscu):
	if type(rscu) == str:
		with open(rscu) as r:
			rscu = json.load(r)
	
	rscu_seq = []
	for i in range(0, len(seq), 3):
		codon = seq[i : i + 3]
		if codon in rscu.keys():
			rscu_seq.append(rscu[codon])
		else:
			rscu_seq.append(None)

	return(rscu_seq)

#def rankSynCodons(rscu):
#	SynCodonRank = SynonymousCodons
#	for aa in SynonymousCodons.keys():
#		codons = SynonymousCodons[aa].keys()

#def translateSeqToCodonRank(seq, rscu):
#	if type(rscu) == str:
#		with open(rscu) as r:
#			rscu = json.load(r)
#
#	rank_seq = []
#	for i in range(0, len(seq), 3):
#		codon = seq[i : i + 3]
#		if codon in SynonymousCodons.keys():


if __name__ == "__main__":
	print(translateSeqToRSCU('ATGGATCCGACTTCTCGACCCGCTATTGTCATCGACAATGGAACTGGGTATACTAAAATGGGATTTGCGGGTAATGTAGAGCCATGTTTTATCTTACCAACAGTAGTTGCAGTAAACGAGTCGTTTTTGAATCAATCCAAGAGTTCTTCTAAGGCAACTTGGCAAACTCAGCACAACGCTGGAGTTGCTGCGGATCTCGATTTCTATATTGGAGATGAAGCTCTTGCGAAATCGCGGTCTAGTAGTACTCATAATCTTCATTATCCAATTGAGCATGGTCAAGTTGAGGATTGGGATGCTATGGAGAGATATTGGCAGCAGTGTATATTTAACTATTTGAGATGTGATCCTGAGGATCATTATTTTCTTCTTACTGAAAGTCCTCTTACTCCTCCTGAAAGTCGAGAATATACTGGAGAAATTTTGTTTGAGACTTTTAATGTTCCTGGGCTTTACATTGCTGTTAATTCAGTTCTTGCACTTGCTGCTGGCTACACAACATCAAAGTGTGAGATGACAGGGGTTGTAGTAGATGTTGGAGATGGGGCAACTCATGTTGTACCTGTTGCAGAAGGTTATGTCATCGGGAGCTGTATCAAATCGATTCCAATTGCTGGCAAAGATGTCACCCTCTTTATCCAGCAACTCATGCGGGAAAGAGGTGAGAATATACCACCAGAAGATTCATTCGATGTAGCCCGTAAAGTGAAGGAAATGTACTGCTACACTTGTTCTGACATCGTAAAGGAGTTTAATAAACATGACAAAGAACCAGCAAAGTATATTAAACAATGGAAAGGTGTAAAGCCAAAGACTGGTGCACCATACACTTGCGACGTGGGATATGAACGATTCCTTGGACCCGAGGTGTTCTTTAATCCAGAGATATACAGCAATGACTTCACAACTACTTTACCAGCTGTGATAGACAAATGTATTCAGTCTGCACCAATTGACACACGAAGAGCTTTATATAAGAATATAGTGTTGTCCGGAGGTTCAACTATGTTCAAAGATTTCGGAAGAAGGTTACAAAGGGACCTCAAGAAGATTGTTGATGCTCGTGTTCTTGCTAATAACGCTCGTACTGGTGGTGAAATTACGTCTCAACCGGTGGAGGTTAATGTCGTGAGCCATCCTGTCCAGAGGTTTGCAGTTTGGTTTGGAGGCTCTGTGCTTTCATCAACTCCTGAGTTTTTCGCGAGTTGCAGAACGAAAGAGGAGTATGAGGAATATGGAGCAAGCATATGCCGCACGAATCCGGTGTTCAAGGGAATGTATTGA', "Ath_actin"))
	print(translateSeqToRSCU('ATGGACGCCGCCTCCCGCCCCGCCGTCGTCATCGACAACGGCACCGGGTACACGAAGATGGGGTTCGCCGGGAACGTGGAGCCCTGCTTCATCACCCCCACCGTCGTCGCCGTCAACGACACCTTCGCCGGCCAGACCAGGGCCAACACCACCAAGGGGAACTGGATGGCGCAGCACAGCGCCGGCGTCATGGCCGATCTCGACTTCTTCATCGGGGAGGACGCCCTGGCCCGCTCCCGCTCCAGCAACACCTACAACCTCAGCTATCCAATTCACAATGGCCAGGTTGAGAATTGGGACACCATGGAGAGGTTCTGGCAGCAGTGCATTTTCAATTACCTGCGGTGTGATCCGGAGGATCACTATTTTCTGCTCACCGAGAGCCCCCTGACTCCTCCCGAGACCCGCGAGTACACCGGGGAGATCATGTTTGAGACTTTCAATGTGCCTGGTCTGTACATAGCGTGCCAGCCGGTTCTTGCCCTCGCAGCAGGATACACCACCACAAAGTGTGAAATGACAGGTGTTGTAGTCGATGTGGGTGATGGGGCTACCCACATTGTTCCTGTTGCTGATGGTTATGTTATAGGAAGCAGCATCAGATCAATTCCAATTACAGGCAAGGATGTTACCCAGTTTATTCAGCAACTCTTAAAGGAAAGAGGTGAGCACATTCCACCAGAAGAATCTTTTGATGTGGCAAGGAGAGTGAAAGAAATGTACTGCTATACATGTTCAGATATTGTGAAGGAATTTAATAAGCATGACAGAGAGCCCAATAAGTACATAAAGCACTGGAGTGGCATCAAACCAAAAACTGGTGCCAAATACACCTGCGACATTGGATATGAACGCTTCCTAGGGCCTGAGATTTTCTTCCACCCTGAGATTTACAACAATGACTTCACCACTCCTTTGCATGTAGTTATTGACAAGTGCATCCAATCATCCCCAATTGACACAAGAAGGGCTCTTTATAAGAATATTGTCTTGTCCGGGGGATCAACTATGTTTAAGGATTTCCACAGGAGGCTGCAGCGGGACCTAAAAAAGATAGTGGATGCACGAGTCCTTGCATCTAATGCTCGACTTGGTGGAGATGCAAAGGCCCAACCTATAGAAGTAAATGTAGTTAGCCATCCTATTCAAAGATATGCTGTCTGGTTTGGTGGTTCTGTGCTTGCCTCTACCGCTGAATTCTACGAGGCTTGCCACACAAAAGCAGAGTATGAAGAGTATGGCGCAAGCATCTGCCGGACGAATCCTGTTTTCAAAGGGATGTACTAA', "Os_actin"))
	#calculateRSCU(sys.argv[1])

