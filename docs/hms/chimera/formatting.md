---
id: formatting
title: Formatting Chimera Output
---

Once you have finished running the Chimera pipeline, you may want to re-format the output. For example, the raw output contains many columns with "None" in the name: those are mostly just for accounting, since they show how many sites were skipped for that feature (for example, some features may be undefined for UNK residues, etc.).

The raw output also does not contain a "Site By Structure" column, which is important for site identification.

The raw output also does not contain any "averages" features, so that needs to be added.

In this repo, the `format_output.py` script accomplishes that task:

```text
julian-ThinkPad-T460:test$ python3 format_output.py --help
Usage: format_output.py [OPTIONS]

Options:
  --features TEXT  Path to the unformatted features file. Required.
  --columns TEXT   A file with the columns to extact/format, one per line
                   (Default: desired_columns.txt)
  --outfile TEXT   The output file (default: chimeraFormat_out_[date].csv)
  --help           Show this message and exit.
```

In the past, the formatting process has been a bit more involved, since I had to carefully match new feature columns with old ones. That past effort is mostly documented in Jupyter notebooks located in the DropBox at `/home/julian/Dropbox (HMS)/Protein_Oxidation/Feature_Engineering/Chimera_Features/Formatting/Scratch`.
