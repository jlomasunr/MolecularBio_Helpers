#!/usr/bin/env python3
import argparse
from Bio import SeqIO

# Currently in progress...

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
    prev_gene_end = 0
    found = False
    # Get GFF entry for gene id
    with open(args.gff) as gff:
        for line in gff:
            if not("#" in line):
                if line.split("\t")[2] == "gene":
                    if id in line:
                        start = int(line.split("\t")[3])-1
                        end = int(line.split("\t")[4])-1
                        chr = line.split("\t")[0]
                        coordinates[id] = {"start":start, "end":end, "chr":chr, "upstream_lim":prev_gene_end}
                        found = True
                    elif found: 
                        coordinates[id]["downstream_lim"] = int(line.split("\t")[3])
                        break
                    else:
                        prev_gene_end = int(line.split("\t")[4])

# Extract sequences from the Genome
for id in coordinates:
    for seq in SeqIO.parse(args.genome, "fasta"):
        if coordinates[id]["chr"] in seq.id:
            if coordinates[id]["start"] - coordinates[id]["upstream_lim"] >= args.p:
                promoter_length = args.p
            else:
                promoter_length = coordinates[id]["start"] - coordinates[id]["upstream_lim"]

            if coordinates[id]["end"] - coordinates[id]["downstream_lim"] >= args.t:
                terminator_length = args.t
            else:
                terminator_length = coordinates[id]["end"] - coordinates[id]["downstream_lim"]

            promoter = seq[coordinates[id]["start"] - promoter_length: coordinates[id]["start"]]
            terminator = seq[coordinates[id]["end"]:coordinates[id]["end"] + terminator_length]
            promoter.id = f"pro{id} {coordinates[id]['chr']}: {coordinates[id]['start'] - promoter_length} - {coordinates[id]['start']}"
            terminator.id = f"ter{id} {coordinates[id]['end']}: {coordinates[id]['start'] - promoter_length} - {coordinates[id]['end'] + terminator_length}"

            with open(args.o, "w") as out:
                SeqIO.write(promoter, out, "fasta")
                SeqIO.write(terminator, out, "fasta")

print(coordinates)