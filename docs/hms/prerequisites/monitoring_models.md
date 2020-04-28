---
id: monitoring_models
title: Maintaining and monitoring modelling jobs
---

This assumes that you already have I-TASSER and/or QUARK jobs running in parallel.

Both of these modelling programs have very high failure rates and they all finish at different times. And, importantly, these jobs often fail without proper error codes, meaning that the job will be marked "COMPLETED" in an `sacct` output, but will not produce models.

So, it is necessary to regularly monitor job status.

## Moving model*.pdb files

The main output from I-TASSER and QUARK jobs are the model*.pdb files (model1.pdb, model2.pdb, model3.pdb). These are the actual protein models. There are other interesting outputs (like B-factor for I-TASSER) that we haven't used for this pipeline, but you might consider using in the future.

So, the best way to check whether a job is done is to look for a model*.pdb file.

To save space, I will often archive a `main` subdirectory after I find a copy model*.pdb files to a seperate location, to save space.

### I-TASSER

model*.pdb files in I-TASSER are generally found directly in the job's subdirectory. So, to search through all potential model*.pdb files in any batch of I-TASSER jobs, I search for `/n/groups/drad/I-TASSER4.3/BatchJobs/*/main/*/model*.pdb`.

After finding one of these model files, I generally rename them with `I-TASSER` (e.g. `I-TASSER_model3.pdb`) and then put them in a folder corresponding to their protein and fragment of origin. In the future, it might be better to name the file with it's header (e.g. `I-TASSER_Q90928_full_protein_model3.pdb`).

### QUARK

model*.pdb files in QUARK are found within the 'cluster' subdirectory. So, to find them, I search for `/n/groups/drad/QUARKmod/BatchJobs/*/main/*/cluster/model*.pdb`.

### General Script

So, we would like to programatically find folders that have finished models, move those models to a centralized location, and then archive those files.

I will move all models to `/n/groups/drad/all_pdb_files`. I would like to archive the jobs in the Silver Lab research file server, located at `/n/files/SysBio/SILVER LAB/`.

```{bash}
#!/usr/bin/env bash
set -e

# Enable globstar for looping through files
shopt -s globstar

# Make an output directory
output_dir="output/archive_run_out_$(date +%Y%m%d_%H%M%S)"
mkdir -p $output_dir

# Glob-based looping is whitespace-safe and recursive
for file in /n/groups/drad/I-TASSER4.3/BatchJobs/*/main/*/model*.pdb; do
        # Put the filename in a log file
        echo $file >>"$output_dir/I-TASSER_models"

        sourceFolder=$(cut -d'/' -f1-9 <<< $file)
        modelNumber=$(cut -d'/' -f10 <<< $file)
        fragment=$(cut -d'/' -f9 <<< $file)
        protein=$(cut -d'_' -f1 <<< $fragment)

        mkdir -p "/home/js741/all_pdb_models/human/${protein}/${fragment}"

        # (rsync -a: archive -c: verify file with checksum)
        rsync -ac $file "/home/js741/all_pdb_models/human/${protein}/${fragment}/I-TASSER_${modelNumber}"

        if [ "$modelNumber" = "model3.pdb" ]
        then
        echo "Found $file : archiving">>"$output_dir/I-TASSER_model3"

        # Compress and save folder in scratch space.
        mkdir -p "/n/scratch2/js741/modeling_archives/I-TASSER/${protein}/${fragment}"
        tar cjPf "/n/scratch2/js741/modeling_archives/I-TASSER/${protein}/${fragment}/archive_12032019.tar.bz" ${sourceFolder} --remove-files
        fi
done

for file in /n/groups/drad/QUARKmod/BatchJobs/*/main/*/cluster/model*.pdb; do
        # Put the filename in a log file
        echo $file >>"$output_dir/QUARK_models"

        sourceFolder=$(cut -d'/' -f1-9 <<< $file)
        modelNumber=$(cut -d'/' -f11 <<< $file)
        fragment=$(cut -d'/' -f9 <<< $file)
        protein=$(cut -d'_' -f1 <<< $fragment)

        mkdir -p "/home/js741/all_pdb_models/human/${protein}/${fragment}"
        rsync -ac $file "/home/js741/all_pdb_models/human/${protein}/${fragment}/QUARK_${modelNumber}"

        if [ "$modelNumber" = "model3.pdb" ]
        then
        echo "Found $file : archiving">>"$output_dir/QUARK_model3"
        # Compress and save folder in scratch space
        mkdir -p "/n/scratch2/js741/modeling_archives/QUARK/${protein}/${fragment}/"
        tar cjPf "/n/scratch2/js741/modeling_archives/QUARK/${protein}/${fragment}/archive_12032019.tar.bz" ${sourceFolder}
        fi
done
```
