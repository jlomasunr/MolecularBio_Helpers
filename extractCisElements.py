#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import math

# Currently in progress...
# Use agat_sp_add_start_and_stop.pl to add start and stop codons to GFF

parser = argparse.ArgumentParser(description="Extract promoters and terminators for Gene ID")
parser.add_argument("-p", help = "Promoter length (default: 1500bp)", default=1500)
parser.add_argument("-t", help="Terminator length (default: 500bp)", default=500)
parser.add_argument("-o", help="Output file name (default: cisElements.fa", default="cisElements.fa")
parser.add_argument("-f", help="Use a text file with gene IDs (default=0)", default=0)
parser.add_argument("id", help="Gene ID or text file of Gene IDs (--f must be 1 to use text file)")
parser.add_argument("genome", help="Genome file (fasta)")
parser.add_argument("gff", help="Genome anntotation file (GFF3)")
args = parser.parse_args()

# Read in ids from text file or use single
if args.f == 1:
    with open(args.id) as file:
        lines = file.readlines()
        ids = [line.rstrip() for line in lines]
else:
    ids = [args.id]

# Get extraction parameters from the GFF
coordinates = {}
for id in ids:
    start = 0
    end = 0
    chr = ''
    prev_gene_end = {"+":0, "-":0}
    found = False
    start_codon = False
    # Get GFF entry for gene id
    with open(args.gff) as gff:
        for line in gff:
            if not("#" in line):
                params = line.split("\t")
                if params[2] == "gene":
                    start = int(params[3])
                    end = int(params[4])
                    chr = params[0]
                    strand = params[6]
                    if id in line:
                        coordinates[id] = {"start":start, "end":end, "chr":chr, "upstream_lim":prev_gene_end[strand], "strand":strand}
                        found = True
                    elif found:
                        if params[6]== coordinates[id]["strand"]:
                            coordinates[id]["downstream_lim"] = start
                            break
                    else:
                        prev_gene_end[strand] = end
                elif found:
                    if params[2] == "start_codon" and not start_codon:
                        coordinates[id]["start_codon"] = int(params[3])
                        start_codon = True
                    elif params[2] == "stop_codon":
                        coordinates[id]["stop_codon"] = int(params[4])

# Extract sequences from the Genome
for id in coordinates:
    # Get promoter length
    if abs(coordinates[id]["start_codon"] - coordinates[id]["upstream_lim"]) >= args.p:
        promoter_length = args.p
    else:
        promoter_length = abs(coordinates[id]["start"] - coordinates[id]["upstream_lim"])
    
    # Get terminator length
    if abs(coordinates[id]["end"] - coordinates[id]["downstream_lim"]) >= args.t:
        terminator_length = args.t
    else:
        terminator_length = abs(coordinates[id]["end"] - coordinates[id]["downstream_lim"])
    
    # Search chromosomes and extract sequence
    for seq in SeqIO.parse(args.genome, "fasta"):
        if coordinates[id]["chr"] in seq.id:
            if coordinates[id]["strand"] == "+":
                promoter = seq[coordinates[id]["start_codon"] - promoter_length: coordinates[id]["start_codon"]-1]
                terminator = seq[coordinates[id]["stop_codon"]+1:coordinates[id]["stop_codon"] + terminator_length]
                promoter.id = f"pro{id}"
                promoter.description = f"| {coordinates[id]['chr']}: {coordinates[id]['start'] - promoter_length} - {coordinates[id]['start']}"
                terminator.id = f"ter{id}" 
                terminator.description=f"| {coordinates[id]['end']}: {coordinates[id]['start'] - promoter_length} - {coordinates[id]['end'] + terminator_length}"
            elif coordinates[id]["strand"] == "-":
                promoter = seq[coordinates[id]["start_codon"]+1: coordinates[id]["start_codon"] + promoter_length].reverse_complement()
                terminator = seq[coordinates[id]["stop_codon"] - terminator_length:coordinates[id]["stop_codon"]-1].reverse_complement()
                promoter.id = f"pro{id} {coordinates[id]['chr']}: {coordinates[id]['start'] - promoter_length} - {coordinates[id]['start']}"
                terminator.id = f"ter{id} {coordinates[id]['end']}: {coordinates[id]['start'] - promoter_length} - {coordinates[id]['end'] + terminator_length}"
            
            with open(args.o, "w") as out:
                SeqIO.write(promoter, out, "fasta")
                SeqIO.write(terminator, out, "fasta")

print(coordinates)