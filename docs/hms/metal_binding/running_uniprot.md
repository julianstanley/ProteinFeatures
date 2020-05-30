---
id: running_uniprot
title: Running uniprotmetal
---

It is pretty straightforward to extract metal binding sites from UniProt.

First, you should grab xml files for each ID of interest from UniProt. To do this, I provided `get_all_uniprot.sh` in this repo:

```bash
#!/usr/bin/env bash
# Usage: bash get_all_uniprot.sh [IDs] [outfolder]
# IDs: A newline-seperated file containing all uniprot IDs to be queried
# outfolder: folder to put the files in! (do not include final)
# Note: downloads in xml format

while read up; do
    wget -O "$2/${up}.xml" "https://www.uniprot.org/uniprot/${up}.xml"
done <$1
```

Then, you can use python to parse those xml files. To do this, I provided `xml_to_binding.py` in this repo.

That script is a bit more involved, usage:

```bash
Usage: python3 xml_to_binding.py [xml_folder]
Where xml_folder is the folder where all of the uniprot XML files are located
```

It includes a few nice features. For example, it converts from UniProt terms to simple metal binding. For example, 'Iron-sulfur (4Fe-4S)' maps to just 'Fe', and "Calcium 2" maps to just "Ca". Here are the full mappings:

```python
METAL_MAPPING = {
    'Iron (heme proximal ligand)': ["Fe"],
    'Zinc 3': ["Zn"],
    'Iron (heme axial ligand)': ["Fe"],
    'Calcium 2; via carbonyl oxygen': ["Ca"],
    'Magnesium': ["Mg"],
    'Calcium 2': ["Ca"],
    'Magnesium 1': ["Mg"],
    'Magnesium or manganese': ["Mg", "Mn"],
    'Calcium; via carbonyl oxygen': ["Ca"],
    'Zinc 2': ["Zn"],
    'Divalent metal cation; catalytic': ["Ca", "Cr", "Co", "Cu",
                                         "Fe", "Pb", "Mn", "Mg", "Ni", "Sn", "Zn"],
    'Iron-sulfur (4Fe-4S)': ["Fe"],
    'Zinc; catalytic': ["Zn"],
    'Calcium 1': ["Ca"],
    'Calcium': ["Ca"],
    'Zinc 4': ["Zn"],
    'Iron; catalytic': ["Fe"],
    'Zinc 1': ["Zn"],
    'Iron (heme distal ligand)': ["Fe"],
    'Calcium 1; catalytic': ["Ca"],
    'Calcium 3; via carbonyl oxygen': ["Ca"],
    'Copper': ["Cu"],
    'Iron-sulfur (4Fe-4S-S-AdoMet)': ["Fe"],
    'Magnesium 2': ["Mg"],
    'Calcium 3': ["Ca"],
    'Zinc': ["Zn"]
}
```

If you run the script as specified, it will print out all uniprot-derived metal binding features.
