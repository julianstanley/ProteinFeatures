---
id: prereq_overview
title: Prerequisites Overview
---

The biggest prereq, especially for people using this package internally, is obtaining a protein structure.

When a structure is not available in the PDB, you need to do two things: (1) decide which fragments of proteins to model, and (2) model those fragments.

To choose which fragments to model, I wrote some python scripts that I'll go over in future sections.

To model fragments, we will use two tools: I-TASSER (for template-directed modelling) and QUARK (for de-novo modelling). Both of these tools are located on O2 at `/n/groups/drad/`.