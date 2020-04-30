#!/usr/bin/env bash
# usage: ./jobrunner_helper.sh [-s scriptdir] [-o outscriptdir] [-j jobtype] [-f filetype]
# example: ./jobrunner_helper.sh -s /n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/scripts_roger_03302020 -o /n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/finished_scripts -j ITASSERjob -f I-TASSER-model3.pdb
# -s scriptdir (REQUIRED) the directory where runscripts are located. You would like to know whether those runscripts ran successfully.
# -o outscriptdir (REQUIRED) where to move scripts that have run successfully.
# -j jobtype (REQUIRED) Some sort of unique identifier located in the JobName of scripts. So, for example, the name of all I-TASSER jobs starts with "ITASSERjob"
# -f filetype (REQUIRED) The name of an output file to search for in /n/groups/drad/all_pdb_models/

# Enable using stars in glob
shopt -s globstar

# Parse args
while getopts s:o:j:f: option
    do
        case "${option}"
        in
        s) SCRIPTDIR=${OPTARG};;
        o) OUTSCRIPTDIR=${OPTARG};;
        j) JOBTYPE=${OPTARG};;
        f) FILETYPE=${OPTARG};;
        esac
    done

# scriptdir is required
if [[ $SCRIPTDIR == "" ]]
    then
    echo "Please specify an input script directory with the -s flag"
    exit 1
fi

# outscriptdir is required
if [[ $OUTSCRIPTDIR == "" ]]
    then
    echo "Please specify an output script directory with the -o flag"
    exit 1
fi

# jobtype is required
if [[ $JOBTYPE == "" ]]
    then
    echo "Please specify a jobtype with the -j flag"
    exit 1
fi

# filetype is required
if [[ $FILETYPE = "" ]]
    then
    echo "Please specify a filetype with the -f flag"
    exit 1
fi

echo "Making $OUTSCRIPTDIR"
mkdir -p $OUTSCRIPTDIR
echo "Getting sacct"
sacct --format "JobName%1000">currently_running_jobs

runcounter=0
movecounter=0

echo "looping through $SCRIPTDIR/*.sh"
for filename in $SCRIPTDIR/*.sh
do
    if [[ ! -e "$filename" ]]; then continue; fi
    echo "Filename: $filename"
    BASENAME="${filename##*/}"
    BASENAME="${BASENAME/.sh/}"
    echo "Basename: $BASENAME"
    exit 0
    if [ $(find /n/groups/drad/all_pdb_models -wholename "*/${BASENAME}/${filetype}" | wc -l) -gt 0 ]
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
