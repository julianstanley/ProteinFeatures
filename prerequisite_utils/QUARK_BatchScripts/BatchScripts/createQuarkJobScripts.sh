#!/usr/bin/env bash
# usage: ./createQuarkJobScripts [-i inputdir] [-o outdir] [-p priority] [-t time] [-n ncores]
# -i inputdir specify an input directory
# -o outdir specify an output directory
# -p priority (default: 'long') specify a slurm priority
# -t time (default: 30-00:00) specify a slurm runtime
# -n ncores (default: 1) specify a number of slurm cores

# Parse the -i input flag
while getopts o:i:p:t:m:n: option
do
case "${option}"
in
i) INPUTDIR=${OPTARG};;
o) OUTDIR=${OPTARG};;
p) PRIORITY=${OPTARG};;
t) TIME=${OPTARG};;
n) NCORES=${OPTARG};;
esac
done


# Check to make sure there's an input directory
if [[ $INPUTDIR == "" ]]
then
	echo $INPUTDIR
	echo "Please specify the input directory with the -i flag"
	exit 0
fi

# Check to make sure there's an output directory
if [[ $OUTDIR == "" ]]
then
	echo $OUTDIR
	echo "Please specify the outptu directory with the -o flag"
	exit 0
fi

# Set default priority
if [[ $PRIORITY == "" ]]
then
	PRIORITY='long'
fi

# Set default time
if [[ $TIME == "" ]]
then
	TIME='20-00:00'
fi


# Set default ncores
if [[ $NCORES == "" ]]
then
	NCORES='1'
fi

# Loop through input directories
for protDir in $INPUTDIR/*/
do
	# Figure out the memory from the sequence length
	length=$(cat ${protDir}/seq.fasta | sed -n 2p | wc -m)
	MEM=$( expr 2000 + 22 \* ${length} )
	MEM="${MEM}M"
	
	# Edit the directory names to work with the script
	protDir=${protDir::-1}
	protName=${protDir/${INPUTDIR}\/''}	
	echo $protName
	SCRIPTFILE="${OUTDIR}/${protName}.sh"

	echo "#!/bin/bash" > $SCRIPTFILE 
	echo "#SBATCH -x compute-e-16-[178-193,198-253]" >>$SCRIPTFILE
	echo "#SBATCH -J QUARKjob_${protName}" >> $SCRIPTFILE
	echo "#SBATCH -p ${PRIORITY}" >> $SCRIPTFILE
	echo "#SBATCH -t ${TIME}" >> $SCRIPTFILE
	echo "#SBATCH -o QUARKjob_${protName}.\%j\.out" >> $SCRIPTFILE
	echo "#SBATCH --mem ${MEM}" >> $SCRIPTFILE
	echo "#SBATCH -n ${NCORES}" >> $SCRIPTFILE
	echo "#SBATCH -e QUARKjob_${protName}.\%j\.err" >> $SCRIPTFILE
	echo "/n/groups/drad/QUARKmod/runQUARK.py -seqname ${protName} -datadir ${protDir} -runstyle serial -homoflag real" >> $SCRIPTFILE

done
