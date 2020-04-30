#!/usr/bin/env bash
set -e

# Enable globstar for looping through files
shopt -s globstar

############################################################################################################################
# Purpose outline/brainstorming
# Overview:
# I want to make it really easy to start and re-start modeling jobs in O2
# The first step in this script will be to see what jobs *have* finished running. 
# Those will have a model3.pdb file. If I see a model3.pdb file, then I know that the job is done. 
# --> In that case, I should (1) copy the pdb file into a standardized location, (2) check to make sure they're identical
# (for 2: use cmp --silent $f1 $f2) (3) archive & move folder with pdb file
# Later: check for all pdb files, delete scripts with completed pdb files
# How to check for running jobs? Right now there are none. Don't we do something with the name, though?
# ^^ Yes, can filter that way
#
############################################################################################################################
output_dir="output/archive_run_out_$(date +%Y%m%d_%H%M%S)"
mkdir -p $output_dir


############################################################################################################################
# Part 1: - Get all .pdb files from folders 
#         - put them in a safe spot in my home directory
#         - compress 
#
#
############################################################################################################################

# Models folder has a format: (Protein) -> (Segment) -> (Models)

# Loop through all model*.pdb files for I-TASSER, copy them to the models folder.
# Looping is whitespace-safe and recursive
for file in /n/groups/drad/I-TASSER4.3/BatchJobs/*/main/*/model*.pdb; do
	# Put the filename in a log file
	echo "Transferring $file"

	sourceFolder=$(cut -d'/' -f1-9 <<< $file)
	modelNumber=$(cut -d'/' -f10 <<< $file)
	fragment=$(cut -d'/' -f9 <<< $file)
	protein=$(cut -d'_' -f1 <<< $fragment)
    
    if [ ! -f "/home/js741/all_pdb_models/human/${protein}/${fragment}/I-TASSER_${modelNumber}" ]; then
        echo "New file: $file"
    fi

	# Make a folder and copy over this pdb structure
	# (rsync -a: archive -z: compress -c: verify file with checksum)
	mkdir -p "/home/js741/all_pdb_models/human/${protein}/${fragment}"
	rsync -az $file "/home/js741/all_pdb_models/human/${protein}/${fragment}/I-TASSER_${modelNumber}"
    
	#rm $file

	#if [ "$modelNumber" = "model3.pdb" ]
	#then
	#echo "Found $file : archiving">>"$output_dir/I-TASSER_model3"
	# Compress and save folder in scratch space
	#mkdir -p "/n/scratch2/js741/modeling_archives/I-TASSER/${protein}/${fragment}"
	#tar cjPf "/n/scratch2/js741/modeling_archives/I-TASSER/${protein}/${fragment}/archive_12032019.tar.bz" ${sourceFolder} --remove-files
	#fi

done

# Loop through all model*.pdb files for QUARK, copy them to the models folder.
# Looping is whitespace-safe and recursive
for file in /n/groups/drad/QUARKmod/BatchJobs/*/*/*/cluster/model*.pdb; do
	# Put the filename in a log file
	echo "Transfering $file (QUARK)"

	sourceFolder=$(cut -d'/' -f1-9 <<< $file)
	modelNumber=$(cut -d'/' -f11 <<< $file)
	fragment=$(cut -d'/' -f9 <<< $file)
	protein=$(cut -d'_' -f1 <<< $fragment)

    if [ ! -f "/home/js741/all_pdb_models/human/${protein}/${fragment}/QUARK_${modelNumber}" ]; then
        echo "New file: $file"
    fi

	# Make a folder and copy over this pdb structure
	# (rsync -a: archive -z: compress -c: verify file with checksum)
	mkdir -p "/home/js741/all_pdb_models/human/${protein}/${fragment}"
	rsync -az $file "/home/js741/all_pdb_models/human/${protein}/${fragment}/QUARK_${modelNumber}"
	#rm $file

	#if [ "$modelNumber" = "model3.pdb" ]
	#then
	#echo "Found $file : archiving">>"$output_dir/QUARK_model3"
	# Compress and save folder in scratch space
	#mkdir -p "/n/scratch2/js741/modeling_archives/QUARK/${protein}/${fragment}/"
	#tar cjPf "/n/scratch2/js741/modeling_archives/QUARK/${protein}/${fragment}/archive_12032019.tar.bz" ${sourceFolder}
	#fi
done

