#!/usr/bin/env python3

'''
Usage: python3 xml_to_binding.py [xml_folder]
Where xml_folder is the folder where all of the uniprot XML files are located
'''

import sys
import glob
from Bio import SeqIO

# Constants
# Residue and then residue number
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
    # Should I do more for divalent metal?
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


def xml_path_to_metals(xml_path):
    ''' Given an xml_path (str) containing Uniprot XML files,
    print out a list of metals for each xml file
    '''
    for xml_file in glob.glob(f"{xml_path}/*.xml"):
        xml_parsed = [x for x in SeqIO.parse(xml_file, "uniprot-xml")][0]
        sequence = xml_parsed.seq
        metal_features = [x for x in xml_parsed.features if x.type ==
                          'metal ion-binding site']
        locations_to_metals = {}
        for metal_feature in metal_features:
            metal_type = metal_feature.qualifiers["description"]
            try:
                metal_types_short = METAL_MAPPING[metal_type]
            except Exception as error:
                print(error)
                raise ValueError(
                    f"No METAL_MAPPING for {metal_type}. Please provide in source code.")

            for metal_type_short in metal_types_short:
                location = metal_feature.location.end
                if location in locations_to_metals:
                    locations_to_metals[location].append(metal_type_short)
                else:
                    locations_to_metals[location] = [metal_type_short]

        for location in locations_to_metals:
            print(
                f'{xml_file.rsplit("/")[-1].rsplit(".xml")[0]}',
                f'_{sequence[location-1]}{location}\t',
                f'{",".join(locations_to_metals[location])}', sep="")


def main():
    ''' Take a xml path from args, print out all metals for each xml file
    '''
    if len(sys.argv) < 2:
        print("Usage: python3 xml_to_binding [path_to_xml_files]")

    else:
        xml_path = sys.argv[1]

    xml_path_to_metals(xml_path)


if __name__ == "__main__":
    main()
