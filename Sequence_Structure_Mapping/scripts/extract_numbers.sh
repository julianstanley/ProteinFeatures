#!/usr/bin/env bash

### Extracts the residue numberings from a given PDB file.
### Usage: ./extract_numbers [pdb-file]

programname=$0

function usage {
    echo "usage: $programname [pdb-file] [output-folder]"
}

if [ $# -eq 0 ]; then
    usage
    exit 1
fi 

name="$1"
name_base="${name##*/}"
numbers_f3=$(grep -P "ATOM\s+[1-9]" $1 | sed 's/\s\s*/ /g' | grep -Po "\s[A-Z][A-Z][A-Z]\s[0-9][0-9]*\s" | cut -f3 -d' ' | uniq | sed ':a;N;$!ba;s/\n/,/g')
numbers_count_f3=$(echo $numbers_f3 | wc -c)

if [[ $numbers_count_f3 -lt 3 ]]; then
    numbers_f4=$(grep -P "ATOM\s+[1-9]" $1 | sed 's/\s\s*/ /g' | grep -Po "\s[A-Z][A-Z][A-Z]\s[A-Z]\s[0-9][0-9]*\s" | cut -f4 -d' ' | uniq | sed ':a;N;$!ba;s/\n/,/g')
    echo $numbers_f4 > "$2/$name_base"
else
    echo $numbers_f3 > "$2/$name_base"
fi
