#!/usr/bin/bash
# Runs the core chimera feature extraction script with all arguments passed
# to this script while also starting the reply log
if [[ $* == *--help* ]]; then
    runHelp="python run_ChimeraFeatureExtraction_core.py --help"
    eval $runHelp
else
    runChimera="chimera --start 'Reply Log' --script 'run_ChimeraFeatureExtraction_core.py $*'"
    echo "Running ${runChimera}"
    eval $runChimera
fi
