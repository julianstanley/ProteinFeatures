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

I will move all models to `/n/groups/drad/all_pdb_files`. I would like to archive the jobs in the Silver Lab research file server, located at `/n/files/SysBio/SILVER LAB/`. However, that directory is only available on the transfer servers (which you can access by ssh'ing into `transfer.rc.hms.harvard.edu` with your eCommonsID). 

So, I archive the run folder (after finding a model3.pdb, since I only collect 3 models) into `/n/groups/drad/archive` and then later move those into a `/n/files/` directory on the transfer server. (Note: I used to use the `scratch` directory as an intermediate for archiving, but we recently increased the storage limit in `/n/groups/drad` so that is no longer necessary.)

The script to archive and save models can be found at `/n/groups/drad/julian/Modelling_Monitoring/transfer_archive_models.sh`, but the format is pretty straightforward, as follows:

```{bash}
#!/usr/bin/env bash
set -e

# Enable globstar for looping through files
shopt -s globstar

# Make an output directories
output_dir="output/archive_run_out_$(date +%Y%m%d_%H%M%S)"
archive_name="archive_$(date +Y%m%d)"
mkdir -p $output_dir

archive_dir="/n/groups/drad/archive"
model_dir="/n/groups/drad/all_pdb_files/human"

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
```

And repeating for QUARK

## Tracking the number of jobs that have completed

Currently (as of April 28), I have a recurring cronjob that checks the status of I-TASSER and QUARK jobs.

The cronjob calls the `/n/groups/drad/all_pdb_models/progress_logs/cronjob.sh` script, which in turn calls a Makefile in the same directory. That makefile first updates all available models (by moving them from `BatchJobs` directories to the `all_pdb_models` directories) and calls a series of relatively simple bash commands to identify completed fragments and models.

For more information, see `/n/groups/drad/all_pdb_models/progress_logs/Makefile`.

To check the number of I-TASSER and QUARK jobs completed over time, I go into `/n/groups/drad/all_pdb_models/progress_logs/updates` and run `head -2 *`. The files in that directory contain a log of the counts of fragments detected.

The `/n/groups/drad/all_pdb_models/progress_logs/old` directory contains more detailed information about those fragments, with folders corresponding to the dates at which those fragments were found. The `total_proteins_count` and `total_fragments_count` are the source tables that I used to update the modelling status.