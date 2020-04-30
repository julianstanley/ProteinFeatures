#!/usr/bin/env bash
# Note: run from the directory containing the .sh runscripts to be tested

# Description: a script used in March and April 2020 to move and run completed and uncompleted jobs.
mkdir -p jobrun_archive
sacct --format "JobName%1000">currently_running_jobs
filetype="I-TASSER_model3.pdb"
jobtype="ITASSERjob"

counter=0
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
    if [ $(cat currently_running_jobs | grep "${jobtype}_${BASENAME}" | wc -l) -gt 0 ]
    then
    echo "Job is already running"
    mv $filename jobs_running_03312020/
    else
    echo "Moving job"
    mv $filename ./scripts_roger_04012020
    ((counter++))
    echo "Counter: $counter"
        if [ $counter -gt 8999 ]; then
        exit 0
        fi
    fi
fi

done
