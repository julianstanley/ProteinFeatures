# Instructions
## Julian Stanley | March 26, 2020

Running DISOPRED is super straightforward. 

1. On O2, locate the DISOPRED program. Right now, it's at `/n/groups/drad/DISOPRED/run_disopred.pl`.
2. Loop through each of your protein FASTA files, and run disopred on that file: e.g. `/n/groups/drad/DISOPRED/run_disopred.pl ${protein} &>${protein}.log`. The output will appear in your current working directory--so, consider `cd`ing into an output directory before running. 
