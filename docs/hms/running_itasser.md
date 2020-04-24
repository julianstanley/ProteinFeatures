---
id: running_itasser
title: Creating and running I-TASSER jobs.
---

All of the tools for running I-TASSER are located in O2 at `/n/groups/drad/I-TASSER4.3/`.

## Running I-TASSER once

Running I-TASSER is fairly straightforward, but it requires a lot of arguments. The main runscript is `/n/groups/drad/I-TASSER4.3/I-TASSERmod/runI-TASSER.pl`.

The perl script takes requires a few arguments to run correctly:

* `-java_home`: Path to a java distribution. We use `/n/app/java/jdk-1.8u112`.
* `-pkgdir`: Path to the main I-TASSER directory. So, `/n/groups/drad/I-TASSER4.3/`.
* `-libdir`: Path to library packages for I-TASSER. So, `/n/groups/drad/I-TASSER4.3/ITLIB`.
* `-seqname`: The sequence name for the protein/fragment you're modelling. For example, `Q9Y5Z9_full_protein`.
* `-datadir`: This is the directory where the data needed to run this I-TASSER job is located. For example, we put the `seq.fasta` for the Q9Y5Z9_full_protein fragment in `/n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/main/Q9Y5Z9_full_protein`.
* `-runstyle`: See the I-TASSER docs for more info. We use `serial`.
* `-homoflag`: See the I-TASSER docs for more info. We use `real`.
* `-nmodel`: How many models to generate. We used `5` and more recently switched to `3`.
* `-LBS`: See the I-TASSER docs for more info. We use `false`.
* `-EC`: See the I-TASSER docs for more info. We use `false`.
* `-GO`: See the I-TASSER docs for more info. We use `false`.
* `-usrnamme`: Your O2 username. For example, Roger's is `rlc18`.

So, here's an example of a full I-TASSER call:

```{bash}
/n/groups/drad/I-TASSER4.3/I-TASSERmod/runI-TASSER.pl -java_home /n/app/java/jdk-1.8u112 -pkgdir /n/groups/drad/I-TASSER4.3 -libdir /n/groups/drad/I-TASSER4.3/ITLIB -seqname Q9Y5Z9_full_protein -datadir /n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/main/Q9Y5Z9_full_protein -runstyle serial -homoflag real -nmodel 3 -LBS false -EC false -GO false -usrname rlc18
```

## Setting up an O2 job for I-TASSER

Although you could run I-TASSER with just the command above (in an interactive job or something), jobs typically take a long time (days to weeks, sometimes over a month). So, you'll probably want to set up a `long` job on O2.

We have been running *all* I-TASSER jobs with `4000M` of memory, `20-00:00` (20 days) of runtime, and a single cpu. You'll see that, for QUARK, we scale our memory requirements based on sequence size. Further troubleshooting is needed to figure automatically figure out how much memory and runtime to allocate to an I-TASSER job.

In addition, in late 2019 the folks at research computing emailed me to tell me that our modelling jobs were overheating some of the server racks. In response, they asked that we only run those jobs (when we are running in high-throughput) on a specific rack of computers.

Putting all of those things together, here is an example runscript for the example protein we used before:

```{bash}
#!/bin/bash
#SBATCH -x compute-e-16-[178-193,198-253]
#SBATCH -J ITASSERjob_Q9Y5Z9_full_protein
#SBATCH -p long
#SBATCH -t 20-00:00
#SBATCH -o /n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/log_files/ITASSERjob_Q9Y5Z9_full_protein.\%j\.out
#SBATCH --mem 4000M
#SBATCH -n 1
#SBATCH -e /n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/log_files/ITASSERjob_Q9Y5Z9_full_protein.\%j\.err
/n/groups/drad/I-TASSER4.3/I-TASSERmod/runI-TASSER.pl -java_home /n/app/java/jdk-1.8u112 -pkgdir /n/groups/drad/I-TASSER4.3 -libdir /n/groups/drad/I-TASSER4.3/ITLIB -seqname Q9Y5Z9_full_protein -datadir /n/groups/drad/I-TASSER4.3/BatchJobs/Final_Push/main/Q9Y5Z9_full_protein -runstyle serial -homoflag real -nmodel 3 -LBS false -EC false -GO false -usrname rlc18
```

## Setting up and running I-TASSER jobs in high-throughput

I created some bash scripts to help you run I-TASSER in high-throughput. You can find them at `/n/groups/drad/I-TASSER4.3/BatchScripts`.

### Create big FASTA file with all sequences to model (csv_to_fasta.sh)

Begin by building a .csv file (I do this on my local computer) that has a line for each fragment you want to model. The first entry in each line should be the name of the fragment (make sure this is unique!) and the second entry should be the sequence of that fragment.

Then, the `/n/groups/drad/I-TASSER4.3/BatchScripts/csv_to_fasta.sh` script is just an `awk` one-liner that you can use to convert that .csv file into FASTA format.

### Split big FASTA file into pieces, making the main input directories for I-TASSER (splitfasta.sh)

Next, you should create a new folder on O2 that will contain all of your runscripts, I make these folders in `/n/groups/drad/I-TASSER4.3/BatchJobs/`. For each folder that I make, I also add a description in the `/n/groups/drad/I-TASSER4.3/BatchJobs/README` file.

So, make a folder within `BatchJobs` with a descriptive name, update the `README`, and then copy your big FASTA-format file into that folder.

Next, you will need to split your big fasta file into a bunch of little ones, one for each job. That's what I made the `/n/groups/drad/I-TASSER4.3/BatchScripts/split_fasta.sh` script for. So, (1) in your BatchJobs subfolder, make a folder that will contain your I-TASSER source directories. I usually call this `main` and then, (2) run the `split_fasta.sh` script, it's usage is like this: `./splitfasta [-i infile] [-o outdir]`.

### From input directories, create the scripts to run individual I-TASSER jobs (createITasserJobScripts.sh)

Next, you need to create runscripts in high-throughput. Create a directory for those scripts--I usually just call it `scripts`.

Then, run the `/n/groups/drad/I-TASSER4.3/BatchScripts/createITasserJobScripts.sh` script, with your `main` directory (the one created in the previous step) as input, and your `scripts` directory as output. The general usage is: `./createITasserJobScripts [-i inputdir] [-o outdir] [-p priority] [-t time] [-n ncores]`.

**TODO note:** This script isn't fully up-to-date yet. There were a couple problems--for example, the input directory was relative, but it should be absolute. And the error and log files should go to a seperate folder to keep everything neat. The `reformat.sh` script shows the different changes that I made post-hoc lately, but haven't gotten around to making in the generator script.

### Run the job scripts (runJobScripts.sh)

Now, you can just run all of the scripts. `find ./scripts -type f -name "*.sh" -exec sbatch {} \;` works fine for this, or you can use my more complicated `runJobScripts.sh` script. The general usage is `./runJobScripts.sh [-s directory with the runscripts]`.

Before you do this, especially if you're running a lot of job, I highly recommend you request an interactive job--just running all the scripts can be taxing.

### Profit

Now you should have lots of I-TASSER jobs in the queue! :)
