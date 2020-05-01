---
id: chimera_overview
title: Chimera Overview
---

Once you have 3D protein structures, you can use UCSF Chimera to extract certain features from them.

It can be a bit involved to extract these features, so I've written a pipeline to do it for you. It requires the UCSF Chimera GUI, which has a built-in python2 interpreter.

Unfortunately, in my hands at least, Chimera (and especially the builtin interpreter) can be very finnicky. So, I included my exact build of UCSF Chimera, built on Ubuntu 19.10, in this repository (`Chimera_Features/UCSF_Chimera64-2019-07-10.zip`). There is also a lab Dell laptop that has the same version of Chimera installed. Different versions can give different results, or sometimes not work alltogether.

## Quick Reference

1. Copy the Chimera_Features folder. Then, create a subdirectory called `pdb` (or something of the sort) and put your source .pdb files in there.

2. Make a subdirectory called `output` or something of the sort.

3. In the root directory, run `bash run_ChimeraFeatureExtraction.sh --pdbdir pdb/ --outdir output/`. This assumes that `chimera` is in your path and that you want to analyze bubble features at 5, 8, 10, and 12A.
