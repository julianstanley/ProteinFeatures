---
id: using_lomets
title: Using LOMETS
---

The LOMETS runscript is in the I-TASSER suite code, specifically at `/n/groups/drad/I-TASSER4.3/I-TASSERmod/runLOMETS.pl`.

To run, LOMETS just requires a path to the I-TASSER package, some library files, java, and the directory where it can find a `.fasta` file containing your sequence of interest.

Specifically, I use these arguments:

* `-pkgdir`: Path to the main I-TASSER directory. So, `/n/groups/drad/I-TASSER4.3/`.
* `-libdir`: Path to library packages for I-TASSER. So, `/n/groups/drad/I-TASSER4.3/ITLIB`.
* `-seqname`: The sequence name (I use uniprot id) for the sequence you're modelling. For example, `Q9Y5Z9`.
* `-java_home`: Path to a java distribution. We use `/n/app/java/jdk-1.8u112`.
* `-runstyle`: See the I-TASSER docs for more info. We use `serial`.
* `-homoflag`: See the I-TASSER docs for more info. We use `real`.
* `-outdir`: Where will LOMETS put intermediate and output files? I usually just use the same directory where I put the `.fasta` file.
* `-usrname`: Your O2 username. For example, Roger's is `rlc18`.

Here's an example full O2 script:

```{bash}
#!/bin/bash
#SBATCH -J LOMETS_O75446
#SBATCH -p short
#SBATCH -t 00-06:00
#SBATCH --mem 1500M
#SBATCH -o /home/js741/CarbonylationSite_Prediction/Human_Proteins/lomets_helpers/LOMETS/main/O75446/logs/lomets_run.out
#SBATCH -e /home/js741/CarbonylationSite_Prediction/Human_Proteins/lomets_helpers/LOMETS/main/O75446/logs/lomets_run.err
/n/groups/drad/I-TASSER4.3/I-TASSERmod/runLOMETS.pl -pkgdir /n/groups/drad/I-TASSER4.3/ -libdir /n/groups/drad/I-TASSER4.3/ITLIB/ -seqname O75446 -datadir /home/js741/CarbonylationSite_Prediction/Human_Proteins/lomets_helpers/LOMETS/main/O75446 -java_home /n/app/java/jdk-1.8u112 -runstyle serial -homoflag real -outdir /home/js741/CarbonylationSite_Prediction/Human_Proteins/lomets_helpers/LOMETS/main/O75446 -usrname js741
```
