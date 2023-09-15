#!/usr/bin/env python

'''
Plasmid comparison for plasmid sequence verification.

Usage: 
ComPlasmid.py [subject.(fa|gbk)] [query.(fa|gbk)] [optional: --transfer]

Purpose: 
Compare a circular DNA sequnce (plasmid) to a reference with potentially
different starting index and orientation. The program re-indexes the query
sequence to match as closely as possible the orietntation of a reference 
plasmid. Sequence features can be transferrred from the reference to the 
query if the reference is provided in genbank format.

Output:
A re-indexed query sequence in fasta and genbank format.

Features are transferred if the subject is provided in 
genbank format and the '--transfer' option is used.
'''

import argparse, re, json, os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, CompoundLocation, SimpleLocation
import parasail

parser = argparse.ArgumentParser(
					prog='ComPlasmid',
					description="Plasmid comparison for plasmid sequence verification.")
parser.add_argument('subject_file', help = "Subject file (fasta or genbank)")
parser.add_argument('query_file', help = "Query file (fasta or genbank)")
parser.add_argument('--transfer', action='store_true', help="Map subject genbank features to query file")
parser.add_argument('--strand', default=None, help="[forward|reverse] Only search one strand for local alignments")
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
#TODO: add option to only search one orientation
if args.strand == "reverse":
	print("Computing reverse local alignments...")
	query_r = query_f.reverse_complement(id=True, description=True, name=True)
	alg_qr = parasail.sg_stats_striped_sat(str(query_r.seq), str(subject.seq), 10, 10, parasail.blosum62)
	query = query_r
	orientation = "reverse"
	loc_alg = alg_qr
elif args.strand == "forward":
	print("Computing forward local alignments...")
	alg_qf = parasail.sg_stats_striped_sat(str(query_f.seq), str(subject.seq), 10, 10, parasail.blosum62)
	query = query_f
	orientation = "forward"
	loc_alg = alg_qf
else:
	print("Computing reverse local alignments...")
	query_r = query_f.reverse_complement(id=True, description=True, name=True)
	alg_qr = parasail.sg_stats_striped_sat(str(query_r.seq), str(subject.seq), 10, 10, parasail.blosum62)
	print("Computing forward local alignments...")
	alg_qf = parasail.sg_stats_striped_sat(str(query_f.seq), str(subject.seq), 10, 10, parasail.blosum62)
	# Identify the best orientation for alignment
	if alg_qf.score > alg_qr.score:
		query = query_f
		orientation = "forward"
		loc_alg = alg_qf
		del query_r
	else:
		query = query_r
		orientation = "reverse"
		loc_alg = alg_qr
		del query_f

# Update query indices to sync the locally aligned region (helps with global alignments for comparison)
dif = int(loc_alg.len_query) - int(loc_alg.end_query)
query.seq = query.seq[loc_alg.end_query:] + query.seq[0:loc_alg.end_query]
query.features = [feat._shift(dif) for feat in query.features]

report["Re-indexing"] = {
	"Orientation" : orientation,
	"Best Alignment": f"End Ref: {loc_alg.end_ref}, End Query: {loc_alg.end_query}",
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
			alg_feat = parasail.sw_stats_striped_sat(str(feature.extract(subject).reverse_complement().seq), str(query.seq), 5, 1, parasail.blosum62)
		else:
			alg_feat = parasail.sw_stats_striped_sat(str(feature.extract(subject).seq), str(query.seq), 5, 1, parasail.blosum62)
		
		matches = alg_feat.matches
		match_pct = round((matches/alg_feat.length)*100, 2)

		quals = feature.qualifiers
		if "note" in quals.keys():
			quals["note"] += [f"Transferred from {os.path.basename(args.subject_file)} by ComPlasmid"]
		else:
			quals["note"] = [f"Transferred from {os.path.basename(args.subject_file)} by ComPlasmid"]
		
		feat_start = alg_feat.end_ref-alg_feat.length
		if feat_start < -1:
			feat_start = alg_feat.len_ref + feat_start
			loc = CompoundLocation([SimpleLocation(feat_start+1, alg_feat.len_ref+1, strand=feature.strand),
								SimpleLocation(0,alg_feat.end_ref+1, strand=feature.strand)])
		else:
			loc = SimpleLocation(feat_start+1, alg_feat.end_ref+1, strand=feature.strand)

		feat = SeqFeature(
			type = feature.type,
			qualifiers = quals,
			location = loc)
		
		query.features += [feat]

		report["Features Transferred"][label] = {
			"Length": len(feature),
			"Identical bases": f"{matches} ({match_pct}%)"
		}


# Ouput re-indexed and annotated query
#TODO: Update locus line with descriptive title
query.annotations["molecule_type"] = "DNA"
query.id = f"{os.path.basename(args.query_file)}.newindex" 
SeqIO.write(query, args.query_file + ".newindex.gbk", "genbank")
SeqIO.write(query, args.query_file + ".newindex.fa", "fasta")
print("Done...\n")
print(f"Re-indexed output files written to: \n\t{args.query_file + '.newindex.gbk'}\n\t{args.query_file + '.newindex.fa'}")

print("############ REPORT ############\n")
print(re.sub("\n\s*\n", "\n", re.sub('[\[\]{},"]', '', json.dumps(report, indent=4))))
