---
id: quick_reference
title: Quick Reference
---

A quick reference for important folders and notes.

## Naming Convention

All headers are named "Protein_FragmentID_LocationStart_LocationEnd", where location is optional.

* `O15347_PF00505_HMG_box_93_161` is a PF00505_HMG_box pfam from the O15347 protein located from residues 93 to 161, inclusive.

* `O15347_full_protein` is the full protein sequence for the O15347 protein.

PDB files are nested in folders with their protein and fragment identifiers.

* The .pdb models for `O15347_PF00505_HMG_box_93_161` are located at `/n/groups/drad/all_pdb_models/human/O15347/O15347_PF00505_HMG_box_93_161`

* The .pdb models for `O15347_full_protein` are located at `/n/groups/drad/all_pdb_models/human/O15347/O15347_full_protein`

## Modelling utility scripts

Utilities to generate pfam-aware fragments to model: `Modelling_Pfam_Selection`, this git repository.

Utilies to run I-TASSER in high-throughput: `/n/groups/drad/I-TASSER4.3/BatchScripts`.

Utilities to run QUARK in high-throughput: `/n/groups/drad/QUARKmod/BatchScripts`.

Utilities to move, archive, and run modelling files: `/n/groups/drad/julian/Modelling_Monitoring`.

Utilites to track the progress of runs: `/n/groups/drad/all_pdb_models/progress_logs`.

## Modelling pipeline

0. Decide what sequences to model, usually using the scripts `Modelling_Pfam_Selection` in the git repository.

1. Create csv file in format "header,sequence"; upload to O2 in a new BatchJobs directory. csv file for I-TASSER should have sequences under 1500aa, csv file for QUARK should be under 800aa.

2. Run `csv_to_fasta.sh`, make a `main` and `scripts` folder within the BatchJobs directory.

3. Run `createITasserJobScripts.sh` and/or `createQuarkJobScripts.sh` with `-i main -o scripts`.

4. Run `runJobScripts.sh` with `-s scripts`.

5. Each day, monitor job progress. Use `/n/groups/drad/julian/Modelling_Monitoring/transfer_archive_models.sh` to centralize all models and archive jobs, and then `/n/groups/drad/all_pdb_models/progress_logs/Makefile` to generate a summary of the completed jobs.

6. TODO: Running new scripts
