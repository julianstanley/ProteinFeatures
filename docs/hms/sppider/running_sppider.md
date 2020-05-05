---
id: running_sppider
title: Running SPPIDER
---

This works best if you are using Outlook as your email service.

1. First, create a rule in your email inbox that will put all SPPIDER output files in the same folder. The subject line of the results emails is always "SPPIDER protein interface recognition results".

2. Then use the `sppider_mech_JAS.pl` script in this repo to run pdb files through the SPPIDER web server. Edit line 32 with your email address, and then run the script on each of your PDB files of interest.

I just use find to loop through all PDB files. So, to submit files in the "pdbs" folder, I would run:

```{bash}
find pdbs/ -type f -exec perl sppider_mech_JAS.pl {} \;
```

It can take some time to submit--maybe around ~10 jobs per minute.

3. Then export emails from the SPPIDER folder into a csv file. You can find that option in the Outlook app for Windows (unfortunately not available on Linux) under File-->Export.

4. Now, you can use the `parseSPPIDER2results.pl` script to print out parsed results. For example:

```{bash}
perl parseSPPIDER2results.pl results_email/SPPIDER_05052020.CSV >parsed_spider_results_05052020.txt
```

The resulting file will have a line for each detected PPI.