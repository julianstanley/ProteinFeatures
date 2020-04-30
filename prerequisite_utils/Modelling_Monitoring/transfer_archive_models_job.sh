#!/bin/bash
#SBATCH -t 5-00:00                         # Runtime of 5 minutes, in D-HH:MM format
#SBATCH -p priority                           # Run in short partition
#SBATCH -o transfer_out_%j.out                 # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=ALL                    # ALL email notification type
#SBATCH --mail-user=julian_stanley@hms.harvard.edu  # Email to which notifications will be sent
/n/groups/drad/julian/Modelling_Monitoring/transfer_archive_models.sh
