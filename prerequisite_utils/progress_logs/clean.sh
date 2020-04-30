#!/usr/bin/env bash
current_date=$(date '+%m_%d_%Y')
outdir="old/$current_date/"
mkdir -p $outdir
find . -maxdepth 1 -type f -name "finalpush_*" -exec mv {} $outdir \;
find . -maxdepth 1 -type f -name "find_*" -exec mv {} $outdir \;
find . -maxdepth 1 -type f -name "pdb*" -exec mv {} $outdir \;
find . -maxdepth 1 -type f -name "update_confirmation" -exec mv {} $outdir \;
find . -maxdepth 1 -type f -name "total*" -exec mv {} $outdir \;
