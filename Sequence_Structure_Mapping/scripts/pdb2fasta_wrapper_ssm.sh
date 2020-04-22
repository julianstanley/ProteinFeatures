#!/usr/bin/env bash

mkdir -p ./structure_fasta

for structure_file in structure_pdb/*.pdb
do
    name_base="${structure_file##*/}"
    name_base="${name_base%.pdb}"
    echo $name_base
    bash ../scripts/pdb2fasta_firstchain_js.sh $structure_file
    mv "${name_base}.fasta" ./structure_fasta
done
