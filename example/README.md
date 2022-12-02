# Contained example

Here is a .zip file containing an example of inputs and outputs for the featurization and alignment code.

Since it's been a couple years since I wrote this and I still don't have detailed documentation, I just wanted to upload this as a reference. 

## Alignments
Example input: `/input/21mers.csv`
Script: python3 ./run_generateAlignments.py --seqfile ./input/21mers.csv --output ./output/alignments.csv
Example output: `/output/alignments/csv`

## Chimera
Example input: `/structures`
Script: bash ./run_ChimeraFeatureExtraction.sh --pdbdir ./structures --radii "5 8 10" --attempts_limit "1"
Example output: `/ChimeraOut_example/individual`
