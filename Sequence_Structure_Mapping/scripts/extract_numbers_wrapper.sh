#!/usr/bin/env bash
mkdir -p structure_numbering
for file in ./structure_pdb/*; do bash ../scripts/extract_numbers.sh $file "structure_numbering"; done
