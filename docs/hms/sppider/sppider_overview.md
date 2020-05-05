---
id: sppider_overview
title: SPPIDER Overview
---

SPPIDER is a web server that can predict sites of protein-protein interaction.

Unfortunately, it is only available as a web server. Forunately, that web server will email you results in plain text.

Just use the `sppider_mech_JAS.pl` script in this repo. Edit line 32 with your email address, and then run the script on each of your PDB files of interest.

I just use find to loop through all PDB files. So, to submit files in the "pdbs" folder, I would run:

```{bash}
find pdbs/ -type f -exec perl sppider_mech_JAS.pl {} \;
```

It can take some time to submit--maybe around ~10 jobs per minute.

