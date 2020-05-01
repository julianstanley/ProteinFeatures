---
id: lomets_overview
title: LOMETS Overview
---

Given a FASTA-format protein sequence, LOMETS finds homologous protein structures that could be used for threading.

LOMETS output files (or equivalent) are required as input for findsitemetal.

## Quick Reference

Helper scripts at: `/n/groups/drad/julian/LOMETS_helpers/`.

Basically, just run `/n/groups/drad/julian/LOMETS_helpers/make_lomets_scripts.sh` in the directory you want your LOMETS jobs to live in.

For example: `bash make_lomets_scripts.sh protein_ids "1500M" "00-06:00" "short"`.
