# Description

Here's a collection of tools to derive a mapping between a protein structure (PDB format file) and a corresponding sequence (FASTA format file).

# Overview

To derive a mapping, you need three things:

(1) Sequences
    (1a) The primary sequence of the PDB-format structure 
    (1b) The primary sequence of the FASTA-format file
(2) PDB Numbering
    (2a) Since PDB files sometimes skip numbers due to gaps/etc., we also
    need a text file with all the site numbers from the PDB-format structre.

After that, a Python script can do the rest. By which I mean it can output a text file
in the following format:

{structure file name}_{residue number from PDB structure}\t{residue number from FASTA sequence}

# Quickstart

```
# Original structure
julian-ThinkPad-T460:04_21_2020$ tree
.
├── sequence_fasta
│   ├── A5A619.fasta
│   ├── P21515.fasta
│   ├── Q46798.fasta
│   ├── Q9RUW8.fasta
│   └── Q9RXJ9.fasta
└── structure_pdb
    ├── A5A619_E10835_deltop10.pdb
    ├── A5A619_E10835_gap30_40.pdb
    ├── A5A619_E10835_mut10_mut50.pdb
    ├── A5A619_E10835.pdb
    ├── P21515_E10907_deltop10.pdb
    ├── P21515_E10907_gap20_30.pdb
    ├── P21515_E10907_mut2_mut82.pdb
    ├── P21515_E10907.pdb
    ├── Q46798_E13810_deltop10.pdb
    ├── Q46798_E13810_gap60_70.pdb
    ├── Q46798_E13810_mut9_mut167.pdb
    ├── Q46798_E13810.pdb
    ├── Q9RUW8_2nvo_A_deltop10.pdb
    ├── Q9RUW8_2nvo_A_gap500_510.pdb
    ├── Q9RUW8_2nvo_A_mut108_mut530.pdb
    ├── Q9RUW8_2nvo_A.pdb
    ├── Q9RXJ9_DR0314_itasser_deltop10.pdb
    ├── Q9RXJ9_DR0314_itasser_gap10_20.pdb
    ├── Q9RXJ9_DR0314_itasser_mut6_mut28.pdb
    └── Q9RXJ9_DR0314_itasser.pdb


#  Will respond with each pdb structure turned into a fasta file, leaving out here for brevity.
julian-ThinkPad-T460:04_21_2020$ bash ../scripts/pdb2fasta_wrapper_ssm.sh 
julian-ThinkPad-T460:04_21_2020$ bash ../scripts/extract_numbers_wrapper.sh 

# Will respond with each structure and its alignment, but leaving out here for brevity
julian-ThinkPad-T460:04_21_2020$ python3 ../scripts/ssm_core.py 

# Final output directory structure
julian-ThinkPad-T460:04_21_2020$ tree
.
├── output_maps
│   ├── individual
│   │   ├── A5A619_E10835.csv
│   │   ├── A5A619_E10835_deltop10.csv
│   │   ├── A5A619_E10835_gap30_40.csv
│   │   ├── A5A619_E10835_mut10_mut50.csv
│   │   ├── P21515_E10907.csv
│   │   ├── P21515_E10907_deltop10.csv
│   │   ├── P21515_E10907_gap20_30.csv
│   │   ├── P21515_E10907_mut2_mut82.csv
│   │   ├── Q46798_E13810.csv
│   │   ├── Q46798_E13810_deltop10.csv
│   │   ├── Q46798_E13810_gap60_70.csv
│   │   ├── Q46798_E13810_mut9_mut167.csv
│   │   ├── Q9RUW8_2nvo_A.csv
│   │   ├── Q9RUW8_2nvo_A_deltop10.csv
│   │   ├── Q9RUW8_2nvo_A_gap500_510.csv
│   │   ├── Q9RUW8_2nvo_A_mut108_mut530.csv
│   │   ├── Q9RXJ9_DR0314_itasser.csv
│   │   ├── Q9RXJ9_DR0314_itasser_deltop10.csv
│   │   ├── Q9RXJ9_DR0314_itasser_gap10_20.csv
│   │   └── Q9RXJ9_DR0314_itasser_mut6_mut28.csv
│   └── map_2020-04-21.csv
├── pickled_intermediates
│   ├── alignments.pickle
│   └── up_to_sequences.pickle
├── sequence_fasta
│   ├── A5A619.fasta
│   ├── P21515.fasta
│   ├── Q46798.fasta
│   ├── Q9RUW8.fasta
│   └── Q9RXJ9.fasta
├── structure_fasta
│   ├── A5A619_E10835_deltop10.fasta
│   ├── A5A619_E10835.fasta
│   ├── A5A619_E10835_gap30_40.fasta
│   ├── A5A619_E10835_mut10_mut50.fasta
│   ├── P21515_E10907_deltop10.fasta
│   ├── P21515_E10907.fasta
│   ├── P21515_E10907_gap20_30.fasta
│   ├── P21515_E10907_mut2_mut82.fasta
│   ├── Q46798_E13810_deltop10.fasta
│   ├── Q46798_E13810.fasta
│   ├── Q46798_E13810_gap60_70.fasta
│   ├── Q46798_E13810_mut9_mut167.fasta
│   ├── Q9RUW8_2nvo_A_deltop10.fasta
│   ├── Q9RUW8_2nvo_A.fasta
│   ├── Q9RUW8_2nvo_A_gap500_510.fasta
│   ├── Q9RUW8_2nvo_A_mut108_mut530.fasta
│   ├── Q9RXJ9_DR0314_itasser_deltop10.fasta
│   ├── Q9RXJ9_DR0314_itasser.fasta
│   ├── Q9RXJ9_DR0314_itasser_gap10_20.fasta
│   └── Q9RXJ9_DR0314_itasser_mut6_mut28.fasta
├── structure_numbering
│   ├── A5A619_E10835_deltop10.pdb
│   ├── A5A619_E10835_gap30_40.pdb
│   ├── A5A619_E10835_mut10_mut50.pdb
│   ├── A5A619_E10835.pdb
│   ├── P21515_E10907_deltop10.pdb
│   ├── P21515_E10907_gap20_30.pdb
│   ├── P21515_E10907_mut2_mut82.pdb
│   ├── P21515_E10907.pdb
│   ├── Q46798_E13810_deltop10.pdb
│   ├── Q46798_E13810_gap60_70.pdb
│   ├── Q46798_E13810_mut9_mut167.pdb
│   ├── Q46798_E13810.pdb
│   ├── Q9RUW8_2nvo_A_deltop10.pdb
│   ├── Q9RUW8_2nvo_A_gap500_510.pdb
│   ├── Q9RUW8_2nvo_A_mut108_mut530.pdb
│   ├── Q9RUW8_2nvo_A.pdb
│   ├── Q9RXJ9_DR0314_itasser_deltop10.pdb
│   ├── Q9RXJ9_DR0314_itasser_gap10_20.pdb
│   ├── Q9RXJ9_DR0314_itasser_mut6_mut28.pdb
│   └── Q9RXJ9_DR0314_itasser.pdb
└── structure_pdb
    ├── A5A619_E10835_deltop10.pdb
    ├── A5A619_E10835_gap30_40.pdb
    ├── A5A619_E10835_mut10_mut50.pdb
    ├── A5A619_E10835.pdb
    ├── P21515_E10907_deltop10.pdb
    ├── P21515_E10907_gap20_30.pdb
    ├── P21515_E10907_mut2_mut82.pdb
    ├── P21515_E10907.pdb
    ├── Q46798_E13810_deltop10.pdb
    ├── Q46798_E13810_gap60_70.pdb
    ├── Q46798_E13810_mut9_mut167.pdb
    ├── Q46798_E13810.pdb
    ├── Q9RUW8_2nvo_A_deltop10.pdb
    ├── Q9RUW8_2nvo_A_gap500_510.pdb
    ├── Q9RUW8_2nvo_A_mut108_mut530.pdb
    ├── Q9RUW8_2nvo_A.pdb
    ├── Q9RXJ9_DR0314_itasser_deltop10.pdb
    ├── Q9RXJ9_DR0314_itasser_gap10_20.pdb
    ├── Q9RXJ9_DR0314_itasser_mut6_mut28.pdb
    └── Q9RXJ9_DR0314_itasser.pdb
```




```

# Creating prerequisites 

1. First, create a folder with today's date and cd into that folder. All of the scripts assume that your working directory is at the same depth as the `scripts` directory. 

This assumes that you start off with (1) a PDB files and (2) a file, `uniprot_ids`, containing each of the canonical uniprot IDs of your PDB files, one ID per lne.

2. Make a folder called `{today's date}/sequence_fasta`. Then, use the `./scripts/get_uniprot.sh` to query uniprot fasta sequences. So, for example, I ran `bash ../../scripts/get_uniprot.sh ../uniprot_ids` while in the `04_20_2020/sequence_fasta` folder. Now, you should have fasta files corresponding to each canonical ID. 

3. Put all of your PDB files in a folder called `{today's date}/structure_pdb`.

Future scripts rely on the names of your files, so make sure to name them with the following scheme:

* For canonical sequence fasta files: {uniprot id}.fasta. So, for example, I have `P00338.fasta` in the `04_20_2020/sequence_fasta` folder. 

* For PDB files: {uniprot id}_{structure name}.pdb. So, for example, I am putting a corresponding `P00338_pdb1i10.pdb` into the `04_17_2020/structure_pdb` folder. 

I created a helper script to help you rename your pdb files, if they are in the format `pdb1i10.pdb` and you want to rename to `P00338_pdb1i10.pdb`. Uniprot
ids are not generally included in PDB files, so you'll need to provide your own comma-delimited mapping between UniProt and PDB structure file. 

^ Remember to back-up your original structures before running this!

## Converting PDB into fasta.

We can use the `./scripts/pdb2fasta_firstchain_js.sh` script to extract a sequence from the the first chain in your PDB file. 

To make this a bit easier on many files, I wrote the `./scripts/pdb2fasta_wrapper_ssm.sh` script.

So, first `cd` into your base directory (the one with the date--in my case, `04_17_2017`).
Then, run `bash ../scripts/pdb2fasta_wrapper_ssm.sh`.

And that should convert all PDB files in `structure_pdb` into fasta files, located in `structure_fasta`. 

## Converting PDB into a numbering file.

We can use the `./scripts/extract_numbers.sh` script to extract sequence numbers from the first chain in your PDB file.

To make this a bit easier on many files, I wrote the `./scripts/extract_numbers_wrapper.sh` script.

Like the last step, you should be in your base directory. Then, run `bash ../scripts/extract_numbers_wrapper.sh`.

And that should convert all PDB files in `structure_pdb` into numbering files, located in `structure_numbering`.

(Note that these numbering files do include all chains, even though we're just using the first chain, so numbering from 
other chains will be ignored)

# Running the main script

Then, you just need to run `python3 ../scripts/ssm_core.py` for the main computation. 

