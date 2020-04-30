#!/usr/bin/env bash

##### Step 1: find all finished pdb files, put paths in a text file
find ../human/ -type f -name *.pdb > find_pdb_all

##### Step 2: find all runscripts in "Final_Push" for I-TASSER and QUARK
find /n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/ -maxdepth 4 -mindepth 2 -type f -name "*.sh" >find_runscripts_finalpush_itasser &

find /n/groups/drad/QUARKmod/BatchJobs/Final_Push/ -maxdepth 4 -mindepth 2 -type f -name "*.sh" >find_runscripts_finalpush_quark &

##### Step 3: Find only model1.pdb lines from the finished pdb files
grep I-TASSER find_pdb_all | grep model1.pdb > find_all_pdb_itasser_model1
grep QUARK find_pdb_all | grep model1.pdb > find_all_pdb_quark_model1


##### Step 4: Get all unique proteins and fragments that have been completed. 

# Cut out unique proteins and fragments from model1 files
cut -d"/" -f4 find_all_pdb_itasser_model1 | sort | uniq >pdb_all_itasser_model1_fragments
cut -d"/" -f3 find_all_pdb_itasser_model1 | sort | uniq >pdb_all_itasser_model1_proteins
cut -d"/" -f4 find_all_pdb_quark_model1 | sort | uniq >pdb_all_quark_model1_fragments
cut -d"/" -f3 find_all_pdb_quark_model1 | sort | uniq >pdb_all_quark_model1_proteins

# Find protein and fragment counts for QUARK and I-TASSER
itasser_fragments=$(wc -l < pdb_all_itasser_model1_fragments)
itasser_proteins=$(wc -l < pdb_all_itasser_model1_proteins)
quark_fragments=$(wc -l < pdb_all_quark_model1_fragments)
quark_proteins=$(wc -l < pdb_all_quark_model1_proteins)

echo "There are ${itasser_fragments} unique fragments from I-TASSER completed, from ${itasser_proteins} proteins."
echo "There are ${quark_fragments} unique fragments from QUARK completed, from ${quark_proteins} proteins."

# Combine counts from QUARK and I-TASSER to get overall counts
cat pdb_all_itasser_model1_fragments pdb_all_quark_model1_fragments | sort | uniq >pdb_all_model1_fragments
cat pdb_all_itasser_model1_proteins pdb_all_quark_model1_proteins | sort | uniq >pdb_all_model1_proteins

total_fragments=$(wc -l < pdb_all_model1_fragments)
total_proteins=$(wc -l < pdb_all_model1_proteins)

echo "There are ${total_fragments} unique fragments in total completed, from ${total_proteins} proteins. Those proteins and fragments are in: 'pdb_all_model1_proteins' and 'pdb_all_model1_fragments', respectively."

##### Step 5: Get all unique proteins and fragments that have been scripts in the works. 

# Repeat the process, this time on scripts in the finalpush files
cut -d"/" -f9 find_runscripts_finalpush_itasser | sed 's/.sh//' | sort | uniq >finalpush_all_itasser_fragments
cut -d"/" -f9 find_runscripts_finalpush_itasser | sed 's/.sh//' | cut -d'_' -f1 | sort | uniq >finalpush_all_itasser_proteins

cut -d"/" -f9 find_runscripts_finalpush_quark | sed 's/.sh//' | sort | uniq >finalpush_all_quark_fragments
cut -d"/" -f9 find_runscripts_finalpush_quark | sed 's/.sh//' | cut -d'_' -f1 | sort | uniq >finalpush_all_quark_proteins

# Find protein and fragment counts for QUARK and I-TASSER, this time with the runscript data
itasser_fragments_scripts=$(wc -l < finalpush_all_itasser_fragments)
itasser_proteins_scripts=$(wc -l < finalpush_all_itasser_proteins)
quark_fragments_scripts=$(wc -l < finalpush_all_quark_fragments)
quark_proteins_scripts=$(wc -l < finalpush_all_quark_proteins)

echo "There are ${itasser_fragments_scripts} unique fragments from I-TASSER with runscripts, from ${itasser_proteins_scripts} proteins."
echo "There are ${quark_fragments_scripts} unique fragments from QUARK with runscripts, from ${quark_proteins_scripts} proteins."

# Combine QUARK and I-TASSER to get overall counts from scripts
cat finalpush_all_quark_fragments finalpush_all_itasser_fragments | sort | uniq >finalpush_all_fragments
cat finalpush_all_quark_proteins finalpush_all_itasser_proteins | sort | uniq >finalpush_all_proteins

total_fragments_scripts=$(wc -l < finalpush_all_fragments)
total_proteins_scripts=$(wc -l < finalpush_all_proteins)

echo "There are ${total_fragments_scripts} unique fragments in total with runscripts, from ${total_proteins_scripts} proteins. Those proteins and fragments are in: 'finalpush_all_proteins' and 'finalpush_all_fragments', respectively."


##### Step 6: Get a combined list of all unique proteins and fragments that are either completed or queued up

cat finalpush_all_fragments pdb_all_model1_fragments | sort | uniq >total_fragments
cat finalpush_all_proteins pdb_all_model1_proteins | sort | uniq >total_proteins

all_fragments_count=$(wc -l <total_fragments)
all_proteins_count=$(wc -l <total_proteins)

echo "There are ${all_fragments_count} unique fragments that have been completed or queued up, from ${all_proteins_count} total proteins. Those proteins and fragments are in: "total_proteins" an "total_fragments", respectively."
