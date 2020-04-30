#!/usr/bin/env bash
# usage: ./createITasserJobScripts [-i inputdir] [-o outdir] [-p priority] [-t time] [-n ncores] 
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

# Set user variable to whoever I am
USER=`whoami`

# Loop through input directories
for protDir in $INPUTDIR/*/
do
	# Figure out the memory from the sequence length
	#length=$(cat ${protDir}/seq.fasta | sed -n 2p | wc -m)
	#MEM=$( expr 1136 + 22 \* ${length} )
	# Keep memory at 4GB until we get a better estimate
	MEM="4000M"
	
	# Edit the directory names to work with the script
	protDir=${protDir::-1}
	# Remove INPUTDIR from protName
	protName=${protDir/${INPUTDIR}\/''}
	# Remove spaces from protName, replace with _
	protName="${protName// /_}"	
	echo $protName
	SCRIPTFILE="${OUTDIR}/${protName}.sh"
	# Remove spaces from SCRIPTFILE, replace with _
	SCRIPTFILE="${SCRIPTFILE// /_}"

	echo "#!/bin/bash" > $SCRIPTFILE 
	echo "#SBATCH -x compute-e-16-[178-193,198-253]" >> $SCRIPTFILE
	echo "#SBATCH -J ITASSERjob_${protName}" >> $SCRIPTFILE
	echo "#SBATCH -p ${PRIORITY}" >> $SCRIPTFILE
	echo "#SBATCH -t ${TIME}" >> $SCRIPTFILE
	echo "#SBATCH -o ITASSERjob_${protName}.\%j\.out" >> $SCRIPTFILE
	echo "#SBATCH --mem ${MEM}" >> $SCRIPTFILE
	echo "#SBATCH -n ${NCORES}" >> $SCRIPTFILE
	echo "#SBATCH -e ITASSERjob_${protName}.\%j\.err" >> $SCRIPTFILE
	echo "/n/groups/drad/I-TASSER4.3/I-TASSERmod/runI-TASSER.pl -java_home /n/app/java/jdk-1.8u112 -pkgdir /n/groups/drad/I-TASSER4.3 -libdir /n/groups/drad/I-TASSER4.3/ITLIB -seqname ${protName} -datadir ${protDir} -runstyle serial -homoflag real -nmodel 3 -LBS false -EC false -GO false -usrname ${USER}" >> $SCRIPTFILE
done
