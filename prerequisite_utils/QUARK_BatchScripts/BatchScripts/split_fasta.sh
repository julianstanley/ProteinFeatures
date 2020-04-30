#!/usr/bin/env bash
# Description: A script to split one large FASTA file into individual sequences
# usage: ./splitfasta [-i infile] [-o outdir]
# -i infile specify an input .fasta file
# -o outdir specify an output directory

# Parse the -o and -i input flags 
while getopts o:i: option
do
case "${option}"
in
o) OUTDIR=${OPTARG};;
i) INFILE=${OPTARG};;
esac
done

# Check to make sure there's an input file
if [[ $INFILE == "" ]]
then
	echo $INFILE
	echo "Please specify an input file with the -i flag"
	exit 0
fi

# Check to make sure that there's an output directory
if [[ $OUTDIR == "" ]]
then
	echo "Please specify an output directory with the -o flag"
	exit 0
fi

# Read the input file
while read line
do
	if [[ ${line:0:1} == '>' ]]
	then
		FILEDIR=${line#>}
		# Remove whitespace, replace with _
		FILEDIR="${FILEDIR// /_}"
		# Remove hypens, replace with _
		FILEDIR="${FILEDIR//-/_}"
		# Remove periods, replace with _
		FILEDIR="${FILEDIR//./_}"

		mkdir -p "${OUTDIR}/${FILEDIR}"
		echo $line > "${OUTDIR}/${FILEDIR}/seq.fasta"
	else
		echo $line >> "${OUTDIR}/${FILEDIR}/seq.fasta"
	fi
done < $INFILE

# Temporary echo output
echo "Input file is $INFILE"
echo "Output directory is $OUTDIR"



