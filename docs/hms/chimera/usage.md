---
id: usage
title: Using Chimera
---

To generate 3D features from Chimera:

```{BASH}
./run_ChimeraFeatureExtraction.sh --pdbdir /home/julian/pdbFiles --radii "5" --attempts_limit '1'
```

## Usage help

```text
Usage: run_generateAlignments.py [OPTIONS]

Options:
  --seqfile TEXT     CSV file with 'sequences' and 'labels' columns
  --matrixfile TEXT  Text file with alignment matrix
  --outfile TEXT     Name of the alignments outfile
  --help             Show this message and exit.
```

```text
Usage: run_ChimeraFeatureExtraction_core.py [OPTIONS]

  Takes in the command line arguments, checks to make sure pdb files were
  passed, and then loops through each pdb file, grabs features from
  chimeraFeatureExtraction.py, formats those features with
  format_single_features, and then writes those features with
  write_all_features() Arguments:     pdbfile (str): Command line argument.
  Path to single pdb file.     pdbdir (str): Command line argument. Path to
  directory of pdb files     outfile (str): Command line argument. Path/name
  of output csv file     radii (str): Space-delineated string of radii at
  which to compute     bubble features. Returns:     Void Effect:     Writes
  a .csv file to "outfile"

Options:
  --pdbfile TEXT         (Optional) PDB file to analyze
  --pdbdir TEXT          (Optional) directory of PDB files
  --outfile TEXT         (Optional) filename to save .csv output file
  --radii TEXT           Space seperated ints representing radii
                         at which to
                         take bubble features
  --attempts_limit TEXT  Number of re-attempts to generate features after
                         failure
  --logfile TEXT         The log file for errors/comments/etc.
  --metals_file TEXT     Path containing the metal binding data file
  --help                 Show this message and exit.
```