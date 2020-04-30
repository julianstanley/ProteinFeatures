#!/usr/bin/env bash
set -e

# Enable globstar for looping through files
shopt -s globstar

# Make an output directories
output_dir="output/archive_run_out_$(date +%Y%m%d_%H%M%S)"
archive_name="archive_$(date +%Y%m%d)"
mkdir -p $output_dir

archive_dir="/n/groups/drad/archive"
model_dir="/n/groups/drad/all_pdb_models/human"

# Glob-based looping is whitespace-safe and recursive
for file in /n/groups/drad/I-TASSER4.3/BatchJobs/*/main/*/model*.pdb; do
        # Put the filename in a log file
        echo $file >>"$output_dir/I-TASSER_models"

        sourceFolder=$(cut -d'/' -f1-9 <<< $file)
        modelNumber=$(cut -d'/' -f10 <<< $file)
        fragment=$(cut -d'/' -f9 <<< $file)
        protein=$(cut -d'_' -f1 <<< $fragment)

        mkdir -p "${model_dir}/${protein}/${fragment}"

        # (rsync -a: archive -c: verify file with checksum)
        rsync -ac $file "${model_dir}/${protein}/${fragment}/I-TASSER_${modelNumber}"

        if [ "$modelNumber" = "model3.pdb" ]
        then
        echo "Found $file : archiving">>"$output_dir/I-TASSER_model3"

        # Compress and save folder in scratch space.
        mkdir -p "${archive_dir}/I-TASSER/${protein}/${fragment}"
        tar cjPf "${archive_dir}/I-TASSER/${protein}/${fragment}/${archive_name}.tar.bz" ${sourceFolder} --remove-files
        fi
done

for file in /n/groups/drad/QUARKmod/BatchJobs/*/main/*/cluster/model*.pdb; do
        # Put the filename in a log file
        echo $file >>"$output_dir/QUARK_models"

        sourceFolder=$(cut -d'/' -f1-9 <<< $file)
        modelNumber=$(cut -d'/' -f11 <<< $file)
        fragment=$(cut -d'/' -f9 <<< $file)
        protein=$(cut -d'_' -f1 <<< $fragment)

        mkdir -p "${model_dir}/${protein}/${fragment}"
        rsync -ac $file "${model_dir}/${protein}/${fragment}/QUARK_${modelNumber}"

        if [ "$modelNumber" = "model3.pdb" ]
        then
        echo "Found $file : archiving">>"$output_dir/QUARK_model3"
        # Compress and save folder in scratch space
        mkdir -p "${archive_dir}/QUARK/${protein}/${fragment}/"
        tar cjPf "${archive_dir}/QUARK/${protein}/${fragment}/${archive_name}.tar.bz" ${sourceFolder}
        fi
done
