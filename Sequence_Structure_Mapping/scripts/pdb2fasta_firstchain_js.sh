#! /bin/bash

# 2012 - P. Poulain

usage="Usage: $0 file.pdb"

#=============================================================================
# input data
#=============================================================================
# check number of arguments
if [ ! $# -eq 1 ]
then
    echo "Argument error" 1>&2
    echo $usage 1>&2
    exit 1
fi

# check first argument is an existing regular file
if [ ! -f $1 ]
then
    echo "$1 is not a regular file" 1>&2
    echo $usage 1>&2
    exit 1
fi

name=$1

#=============================================================================
# functions
#=============================================================================
# list chains in PDB
list_chain() {
    awk '/^ATOM/ && $3 == "CA" {print $5}' | uniq
}

# extract residue sequence from ATOM lines
# take residue from CA atom since
# there is one CA atom per residue
extract_seq_chain() {
    awk -v ch=$1 '/^ATOM/ && $3 == "CA" && $5 == ch {print $4}'
}

# convert newline by space
# for sed lowers this works too: sed ':a;N;$!ba;s/\n/ /g'
remove_newline() {
    tr '\n' ' ' 
}

# convert 3-letter residue code to 1-letter
convert_aa() {
    sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g;s/UNK/X/g;s/ASX/B/g;s/Glx/Z/g;s/Xle/J/g;s/Sec/U/g;s/Pyl/O/g'
}

# remove space between residues
remove_space() {
    sed 's/ //g' 
}

# split fasta sequence at 60 characters (easier to read than 80)
split_60() {
    fold -w 60
}



#=============================================================================
# list chain in PDB file
#=============================================================================
chains=$(cat $name | list_chain)

#=============================================================================
# try to extract a sequence
#=============================================================================
chain=(${chains[0]})
    sequence=$(cat $name | extract_seq_chain  $chain | remove_newline | convert_aa | remove_space)
    sequence_length=$(echo $sequence | wc -l)
    
    if [[ $sequence_length -lt 5 ]]
    then
        sequence=$(cat $name | awk '/^ATOM/ && $3 == "CA" {print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g;s/UNK/X/g;s/ASX/B/g;s/Glx/Z/g;s/Xle/J/g;s/Sec/U/g;s/Pyl/O/g' | sed 's/ //g')
    fi

    if [[ -n $sequence ]]
    then
        name_base="${name##*/}"
        name_base="${name_base%.pdb}"
        size=$(echo $sequence | wc -c)
        size=$((size-1))
        echo ">${name%.pdb} | chain $chain | $size aa" >>"${name_base}.fasta"
        echo $sequence | split_60 >> "${name_base}.fasta"
    fi

