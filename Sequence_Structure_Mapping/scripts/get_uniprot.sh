#!/usr/bin/env bash
# usage: bash get_uniprot.sh [file-name]
# file-name: A path pointing to a file that has one uniprot ID per line
# For each ID, this script will download FASTA sequences from uniprot 
# corresponding to that ID to your current directory.

while read up; do wget -O ${up}.fasta https://www.uniprot.org/uniprot/$up.fasta; done<$1
