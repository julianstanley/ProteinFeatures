#!/usr/bin/env python3
# This script takes in five inputs: (
# (1) SPPIDER output from parseSPPIDER2results.pl, 
# (2) A directory for protein numberings, 
# (3) a directory for protein primary sequences,
# (4) a directory containing the structures to analyze, 
# (5) an output file name

# Here's a hastily written (in a few minutes on my last day) script to convert the output
# of parseSPPIDER2results.pl into the standard format that the Chimera pipeline takes for
# SPPIDER results. It has a few limitations--mostly that it requires that there be 
# files with the FASTA sequence of a structure file and files with the associated indicies
# of that FASTA sequence. I have scripts that can make those files, but imperfectly--
# for example, some PDB formats are very weird and break from those scripts. 
# This parser will warn about those cases, but not break.

import sys
import glob

sppider_output = sys.argv[1]
numberings = sys.argv[2]
sequences = sys.argv[3]
structures = sys.argv[4]
out_filename = sys.argv[5]

# Build a dictionary of SPPIDER hits
hits = {}
with open(sppider_output) as spout:
    for line in spout:
        model, loc = line.strip().rsplit("_", 1)
        aa = loc[0]
        number = loc[1:]
        if model in hits:
            hits[model][number] = aa
        else: 
            hits[model] = {}

with open(out_filename, "w") as outfh:
    success = 0
    fail = 0
    for filename in glob.glob(f'{structures}/*.pdb'):
        basename = filename.rsplit("/", 1)[-1].replace("//", "/")
        numbering_filename = f'{numberings}/{basename}'.replace("//", "/")
        sequence_filename = f'{sequences}/{basename}'.replace("//", "/").replace(".pdb", ".fasta")
        
        with open(numbering_filename) as numbering_fh:
            numbering = numbering_fh.read().strip().split(",")

        with open(sequence_filename) as sequence_fh:
            sequence = "".join([x.strip() for x in sequence_fh.readlines()[1:]])
            
        if len(numbering) != len(sequence):
            print(f"Warning: cannot properly parse: {basename}, skipping.")
            fail += 1
            continue

        if basename in hits:
            for i in range(len(numbering)):
                number = numbering[i]
                aa = sequence[i]
                if number in hits[basename]:
                    if hits[basename][number] == aa:
                        outfh.write(f"{basename},{aa},{number},1\n")
                    else:
                        print(f"Warning: unable to properly match: {basename}, assuming no hits")
                        fail += 1
                        for j in range(i, len(numbering)):
                            number = numbering[j]
                            aa = sequence[j]
                            outfh.write(f"{basename},{aa},{number},0\n")
                        break
            success += 1
        else:
            for i in range(len(numbering)):
                number = numbering[i]
                aa = sequence[i]
                outfh.write(f"{basename},{aa},{number},0\n")
            success += 1

    print(f"Successfully parsed {success} proteins, but failed to parse {fail}.")




