#!/usr/bin/env bash
# Usage: bash get_all_uniprot.sh [IDs] [outfolder]
# IDs: A newline-seperated file containing all uniprot IDs to be queried
# outfolder: folder to put the files in! (do not include final)
# Note: downloads in xml format

while read up; do
    wget -O "$2/${up}.xml" "https://www.uniprot.org/uniprot/${up}.xml"
done <$1
