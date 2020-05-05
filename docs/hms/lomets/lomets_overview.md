---
id: lomets_overview
title: LOMETS Overview
---

Given a FASTA-format protein sequence, LOMETS finds homologous protein structures that could be used for threading.

LOMETS output files (or equivalent) are required as input for findsitemetal.

## Quick Reference

Helper scripts at: `/n/groups/drad/julian/LOMETS_helpers/`.

Basically, just run `/n/groups/drad/julian/LOMETS_helpers/make_lomets_scripts_full.sh` in the directory you want your LOMETS jobs to live in and give a list of uniprot IDs as input.

For example: `bash make_lomets_scripts.sh protein_ids "1500M" "00-06:00" "short"`.

That script will download all uniprot IDs from UniProt directly. If you suspect that some of your IDs can't be downloaded from uniprot, then:

* Use the `/n/groups/drad/julian/LOMETS_helpers/get_ups_from_uniprot.sh` to download all IDs that _can_ be found on uniprot
* Write a custom script to make job folders for the other sequences
* Make a list of all IDs that you want to make scripts for
* Run `/n/groups/drad/julian/LOMETS_helpers/make_lomets_scripts.sh` in the same manner.
