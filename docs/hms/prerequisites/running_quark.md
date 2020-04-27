---
id: running_quark
title: Creating and running QUARK jobs.
---

All of the tools for running QUARK are located in O2 at `/n/groups/drad/QUARKmod`.

## Summary/Quickstart for most cases

1. (On O2) create a folder for your set of jobs in `/n/groups/drad/QUARKmod/BatchJobs`. Add a job description to the README.

2. (On local computer) Create csv file with all fasta headers and sequences, scp to O2.

3. (This step and on are all on O2) Run `bash /n/groups/drad/QUARKmod/BatchScripts/csv_to_fasta.sh your_sequences.csv`

4. In your job folder (for example, `/n/groups/drad/QUARKmod/BatchJobs/myjobs`), run`mkdir main` and `mkdir scripts`.

5. Run `bash/n/groups/drad/QUARKmod/BatchScripts/createQuarkJobScripts.sh -i main -o scripts`.

6. Run `bash /n/groups/drad/QUARKmod/BatchScripts/runJobScripts.sh -s scripts`.

## Running QUARK once

Running QUARK is fairly straightforward. The main runscript is `/n/groups/drad/QUARKmod/runQUARK.py`.

The python script requires a few arguments to run correctly:

* `-seqname`: The sequence name for the protein/fragment you're modelling. For example, `Q9Y5Z9_full_protein`.
* `-datadir`: This is the directory where the data needed to run this QUARK job is located. For example, we put the `seq.fasta` for the Q9Y5Z9_full_protein fragment in `/n/groups/drad/QUARKmod/BatchJobs/Final_Push/main/Q9Y5Z9_full_protein`.
* `-runstyle`: See the QUARK docs for more info. We use `serial`.
* `-homoflag`: See the QUARK docs for more info. We use `real`.

So, here's an example of a full QUARK call:

```{bash}
/n/groups/drad/QUARKmod/runQUARK.py -seqname Q9Y5Z9_full_protein -datadir main/Q9Y5Z9_full_protein -runstyle serial -homoflag real
```

## Setting up an O2 job for QUARK

Although you could run QUARK with just the command above (in an interactive job or something), jobs typically take a long time (days to weeks, sometimes over a month). So, you'll probably want to set up a `long` job on O2.

We have been running all QUARK jobs with `20-00:00` (20 days) of runtime and a single cpu. However, memory is dynamically allocated based on sequence size: equal to `2000M + 22 * length`, where `length` is the sequence length, in AAs. That equation was determined empirically by regressing the used memory of a few thousand previous QUARK runs. A similar analysis could be done with runtime, too (if you have time, you should do this and update these docs!).

In addition, in late 2019 the folks at research computing emailed me to tell me that our modelling jobs were overheating some of the server racks. In response, they asked that we only run those jobs (when we are running in high-throughput) on a specific rack of computers.

Putting all of those things together, here is an example runscript for the example protein we used before:

```{bash}
#!/bin/bash
#SBATCH -x compute-e-16-[178-193,198-253]
#SBATCH -J QUARKjob_Q9Y5Z9_full_protein
#SBATCH -p long
#SBATCH -t 20-00:00
#SBATCH -o QUARKjob_Q9Y5Z9_full_protein.\%j\.out
#SBATCH --mem 9480M
#SBATCH -n 1
#SBATCH -e QUARKjob_Q9Y5Z9_full_protein.\%j\.err
/n/groups/drad/QUARKmod/runQUARK.py -seqname Q9Y5Z9_full_protein -datadir main/Q9Y5Z9_full_protein -runstyle serial -homoflag real
```

## Setting up and running QUARK jobs in high-throughput

I created some bash scripts to help you run QUARK in high-throughput. You can find them at `/n/groups/drad/QUARKmod/BatchScripts`.

### Create big FASTA file with all sequences to model (csv_to_fasta.sh)

Begin by building a .csv file (I do this on my local computer) that has a line for each fragment you want to model. The first entry in each line should be the name of the fragment (make sure this is unique!) and the second entry should be the sequence of that fragment.

Then, the `/n/groups/drad/QUARKmod/BatchScripts/csv_to_fasta.sh` script is just an `awk` one-liner that you can use to convert that .csv file into FASTA format.

### Split big FASTA file into pieces, making the main input directories for QUARK (splitfasta.sh)

Next, you should create a new folder on O2 that will contain all of your runscripts, I make these folders in `/n/groups/drad/QUARKmod/BatchJobs/`. For each folder that I make, I also add a description in the `/n/groups/drad/QUARKmod/BatchJobs/README` file.

So, make a folder within `BatchJobs` with a descriptive name, update the `README`, and then copy your big FASTA-format file into that folder.

Next, you will need to split your big fasta file into a bunch of little ones, one for each job. That's what I made the `/n/groups/drad/QUARKmod/BatchScripts/split_fasta.sh` script for. So, (1) in your BatchJobs subfolder, make a folder that will contain your QUARK source directories. I usually call this `main` and then, (2) run the `split_fasta.sh` script, it's usage is like this: `./splitfasta [-i infile] [-o outdir]`.

### From input directories, create the scripts to run individual QUARK jobs (createITasserJobScripts.sh)

Next, you need to create runscripts in high-throughput. Create a directory for those scripts--I usually just call it `scripts`.

Then, run the `/n/groups/drad/QUARKmod/BatchScripts/createQuarkJobScripts.sh` script, with your `main` directory (the one created in the previous step) as input, and your `scripts` directory as output. The general usage is: `./createITasserJobScripts [-i inputdir] [-o outdir] [-p priority] [-t time] [-n ncores]`.

**TODO note:** This script isn't fully up-to-date yet. There were a couple problems--for example, the input directory was relative, but it should be absolute. And the error and log files should go to a seperate folder to keep everything neat.

### Run the job scripts (runJobScripts.sh)

Now, you can just run all of the scripts. `find ./scripts -type f -name "*.sh" -exec sbatch {} \;` works fine for this, or you can use my more complicated `runJobScripts.sh` script. The general usage is `./runJobScripts.sh [-s directory with the runscripts]`.

Before you do this, especially if you're running a lot of job, I highly recommend you request an interactive job--just running all the scripts can be taxing.

### Profit

Now you should have lots of QUARK jobs in the queue! :)
