---
id: ssm_overview
title: Sequence-Structure Mapping Overview
---

To properly use this pipeline, you need to know the mapping between the canonical protein sequence and the protein structure that you are using.

For example, Lysine #4 in the UniProt sequence may not necessarily be at position #4 in a corresponding protein structure.

In the past, we have annotated those mappings by-hand, which is very time consuming.

In this repo, I include scripts to automate that process by aligning a "canonical" uniprot sequence with a sequence derived from a PDB-format structure file.