---
id: chimera_overview
title: Chimera Overview
---

Once you have 3D protein structures, you can use UCSF Chimera to extract certain features from them.

It can be a bit involved to extract these features, so I've written a pipeline to do it for you. It requires the UCSF Chimera GUI, which has a built-in python2 interpreter.

## Versions

Unfortunately, in my hands at least, Chimera (and especially the builtin interpreter) can be very finnicky. So, I included my exact build of UCSF Chimera, built on Ubuntu 19.10, on the HMS Dropbox (in `Protein_Oxidation/Feature_Engineering/Chimera_Features/Chimera_Backup`). There is also a lab Dell laptop that has the same version of Chimera installed. Different versions can give different results, or sometimes not work alltogether.

**Chimera version: alpha version 1.14 (build 42046) 2019-07-10 00:58:14 UTC Platform: Linux64. Windowing system: x11.**

## Quick Reference

1. Copy the Chimera_Features folder. Then, create a subdirectory called `pdb` (or something of the sort) and put your source .pdb files in there.

2. Make a directory called `input` that will contain source text files for DISOPRED and SPPIDER.

3. Make a subdirectory called `output` or something of the sort.

4. In the root directory, run `bash run_ChimeraFeatureExtraction.sh --pdbdir pdb/ --outdir output/`. This assumes that `chimera` is in your path and that you want to analyze bubble features at 5, 8, 10, and 12A.

## Imports

### Abbreviated List

* chimera
* click
* collections
* datetime
* os
* pandas
* sys
* time
* traceback
* re
* math

Additional dependencies for pandas and click:

Click==7.0
pandas==0.25.0

* numpy [required: >=1.13.3, installed: 1.16.2]
* python-dateutil [required: >=2.6.1, installed: 2.7.3]
* pytz [required: >=2017.2, installed: 2019.2]

### Verbose List

```{text}
run_ChimeraFeatureExtraction_core:
  chimera.runCommand
  click
  collections
  datetime
  os
  pandas
  src.chimeraFeatureExtraction
  sys
  time
  traceback
src.__init__:
  
src.chimeraFeatureExtraction:
  chimera.runCommand
  chimera.selection
  re
  src.ExtendedChimera.ChimeraFeatures
  src.ExtendedChimera.ChimeraUtils
  src.ExtendedChimera.ExtendedAtom
  src.ExtendedChimera.MetalAtom
  src.ExtendedChimera.Utils
src.ExtendedChimera.MetalAtom:
  
src.ExtendedChimera.ExtendedAtom:
  chimera
  math
  re
  src.ExtendedChimera.ChimeraUtils
  src.ExtendedChimera.Utils
src.ExtendedChimera.__init__:
  
src.ExtendedChimera.Utils:
  pandas
src.ExtendedChimera.ChimeraFeatures:
  chimera.runCommand
  chimera.selection
  math
  re
  src.ExtendedChimera.ChimeraUtils
  src.ExtendedChimera.ExtendedAtom
  src.ExtendedChimera.Utils
src.ExtendedChimera.ChimeraUtils:
  AddH
  DockPrep.prefs.INCOMPLETE_SC
  DockPrep.prefs.defaults
  DockPrep.prep
  chimera
  chimera.dialogs
  chimera.runCommand
  src.ExtendedChimera.Utils
```
