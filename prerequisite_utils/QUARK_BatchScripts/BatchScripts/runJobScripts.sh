#!/usr/bin/env bash
# usage: TODO
# Parse the -s script flag
while getopts s: option
do
case "${option}"
in
s) SCRIPTDIR=${OPTARG};;
esac
done

# Check to make sure there's a script directory
if [[ $SCRIPTDIR == "" ]]
then
	echo $SCRIPTDIR
	echo "Please specify the directory containing script files with the -s flag"
	exit 0
fi


for file in ${SCRIPTDIR}/*
do
	sbatch $file
done
