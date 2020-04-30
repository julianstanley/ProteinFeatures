#!/usr/bin/env bash
# Note: run from the directory containing the .sh runscripts to be tested

mkdir -p jobrun_archive
sacct --format "JobName%1000">currently_running_jobs
filetype="I-TASSER_model3.pdb"
jobtype="ITASSERjob"

runcounter=0
movecounter=0
for filename in scripts/*.sh
do
BASENAME=$(echo $filename | sed 's/.sh//' | sed 's/scripts\///')
echo $BASENAME
if [ $(find /home/js741/all_pdb_models -wholename "*/${BASENAME}/${filetype}" | wc -l) -gt 0 ]
then
    echo "Found, moving"
    mv $filename jobrun_archive/
else
    echo "Not Found, running script, if applicable"
    if [ $(grep "${jobtype}_${BASENAME}" currently_running_jobs | wc -l) -gt 0 ]
    then
    echo "Job is already running, moving to 'running' folder"
    mv $filename jobs_running_03312020/
    ((movecounter++))
    echo "Movecounter: ${movecounter}"
    else
    echo "Running job"
    sbatch $filename
    ((runcounter++))
    echo "Runcounter: ${runcounter}" 
    fi
fi

done
