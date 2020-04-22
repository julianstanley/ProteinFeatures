#!/usr/bin/env bash
# usage: bash rename_pdb [mapping-file]
# mapping-file should be comma-delim. uniprot,structure
# Should be run in the same directory as the structure files.
# renames files inplace

while read file
do
    uniprot=$(echo $file | cut -f1 -d',')
    structure_file=$(echo $file | cut -f2 -d',')
    name_base="${structure_file##*/}"
    name_base="${name_base%.pdb}"
    mv ./$structure_file ./${uniprot}_${structure_file}
done<$1

