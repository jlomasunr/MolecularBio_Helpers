#!/usr/bin/env python

'''
Plasmid comparison for plasmid sequence verification.\n

Usage: 
ComPlasmid.py [subject.(fa|gbk)] [query.(fa|gbk)] [optional: --transfer]

Purpose: 
To compare a circular DNA sequnce (plasmid) to a reference with potentially
different starting index and orientation. The query sequence is re-indexed 
to match as closely as possible the orietntation of a reference plasmid. 
Sequence features can be transferrred from the reference to the query if
the reference is provided in genbank format.

Output:
A re-indexed query sequence in fasta and genbank format.

Features are transferred if the subject is provided in 
genbank format and the '--transfer' option is used.
'''

import argparse, re, json
from Bio import SeqIO
from Bio import pairwise2
from Bio import SeqFeature

parser = argparse.ArgumentParser(
					prog='ComPlasmid',
					description="Plasmid comparison for plasmid sequence verification.")
parser.add_argument('subject_file', help = "Subject file (fasta or genbank)")
parser.add_argument('query_file', help = "Query file (fasta or genbank)")
parser.add_argument('--transfer', action='store_true', help="Map subject features to query file")
args = parser.parse_args()

def parseFileType(x:str):
	if re.search("\.gbk$", x):
		return "genbank"
	else:
		return "fasta"

subject_format = parseFileType(args.subject_file)
query_format = parseFileType(args.query_file)

# Read subject and query seqs into memory
print("\nReading input files...")
subject = SeqIO.to_dict(SeqIO.parse(args.subject_file, subject_format))
subject = subject[list(subject)[0]]

query_f = SeqIO.to_dict(SeqIO.parse(args.query_file, query_format))
query_f = query_f[list(query_f)[0]]
query_r = query_f.reverse_complement(id=True, description=True, name=True)

report = {
	"Query" : {
		"file": args.query_file,
		"format": query_format,
		"length": f"{len(query_f.seq)} bp"
		},
	"Reference": {
		"file": args.subject_file,
		"format": subject_format,
		"length": f"{len(subject.seq)} bp"
		}
	}

# Compute local alignments with large mismatch, gap opening, and gap extension 
# costs for both sequence orientations
print("Computing local alignments...")
alg_qr = pairwise2.align.localms(subject.seq, query_r.seq, 5, -100, -100, -50)[0]
alg_qf = pairwise2.align.localms(subject.seq, query_f.seq, 5, -100, -100, -50)[0]

# Identify the best orientation for alignment
if alg_qf.score > alg_qr.score:
	starts = re.findall("\d+", pairwise2.format_alignment(*alg_qf))[:2]
	query = query_f
	orientation = "forward"
	loc_alg = alg_qf
	del query_r
else:
	starts = re.findall("\d+", pairwise2.format_alignment(*alg_qr))[:2]
	query = query_r
	orientation = "reverse"
	loc_alg = alg_qr
	del query_f

# Update query indices to sync the locally aligned region (helps with global alignments for comparison)
dif = int(starts[1]) - int(starts[0])
query.seq = query.seq[dif:] + query.seq[0:dif]
query.features = [feat._shift(-1*dif) for feat in query.features]

report["Re-indexing"] = {
	"Orientation" : orientation,
	"Best Alignment": f"Start: {loc_alg.start}, End: {loc_alg.start}",
	"Shift": dif
	}

# Map subject features to query via local alignments
# TODO: parallelize https://superfastpython.com/multiprocessing-pool-for-loop/
print("Transferring features...")
report["Features Transferred"] = {}
if args.transfer & (subject_format == "genbank"):
	for feature in subject.features:
		label = feature.qualifiers["label"][0] if "label" in feature.qualifiers.keys() else "Misc. Feature"
		
		if feature.strand == -1:
			alg_feat = pairwise2.align.localms(feature.extract(subject).reverse_complement().seq, query.seq, 1, -1, -1, -1)[0]
		else:
			alg_feat = pairwise2.align.localms(feature.extract(subject).seq, query.seq, 1, -1, -1, -1)[0]
		
		matches = pairwise2.format_alignment(*alg_feat).count("|")
		match_pct = round((matches/len(feature))*100, 2)

		quals = feature.qualifiers
		if "note" in quals.keys():
			quals["note"] += [f"Transferred from {args.subject_file} by ComPlasmid"]
		else:
			quals["note"] = [f"Transferred from {args.subject_file} by ComPlasmid"]

		feat = SeqFeature.SeqFeature(
			type = feature.type,
			qualifiers = quals,
			location = SeqFeature.FeatureLocation(start = alg_feat.start, end = alg_feat.end)
			)
		
		query.features += [feat]

		report["Features Transferred"][label] = {
			"Length": len(feature),
			"Identical bases": f"{matches} ({match_pct}%)"
		}


# Ouput re-indexed and annotated query
query.annotations["molecule_type"] = "DNA"
SeqIO.write(query, args.query_file + ".newindex.gbk", "genbank")
SeqIO.write(query, args.query_file + ".newindex.fa", "fasta")
print("Done...\n")
print(f"Re-indexed output files written to: \n\t{args.query_file + '.newindex.gbk'}\n\t{args.query_file + '.newindex.fa'}")

print("############ REPORT ############\n")
print(re.sub("\n\s*\n", "\n", re.sub('[\[\]{},"]', '', json.dumps(report, indent=4))))