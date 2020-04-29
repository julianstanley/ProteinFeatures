import pytest
from src import protein_modeling
from src.protein_modeling import *
import pickle
import copy

# TODO: You need a lot more tests up in here, young padawan. You should do that when you have time

# Setup data
import os
print(os.getcwd())
with open("data/pickle/pickle_07232019_pfam_annotation.txt", "rb") as fp:   # Unpickling
    pfam_annotation = pickle.load(fp)

with open("data/pickle/pickle_07232019_crystal_annotation.txt", "rb") as fp:
    crystal_annotation = pickle.load(fp)

with open("data/pickle/pickle_07232019_protein_crystal.txt", "rb") as fp:
    protein_crystal = pickle.load(fp)

lookup_uniprot_pfams = list_to_lookup(pfam_annotation, "UniProtKB")
lookup_uniprot_crystal = list_to_lookup(protein_crystal, "Protein")
lookup_uniprot_crystals = list_to_lookup(crystal_annotation, "UniProtKB")

protein = 'P16403'

uniprot_crystal_protein = lookup_uniprot_crystal[protein]
full_sequence = uniprot_crystal_protein[0]["Sequence"]
if uniprot_crystal_protein[0]["Crystal"] == "No":
    crystal_dicts = []
else:
 crystal_dicts = lookup_uniprot_crystals[protein]
    
try:
 pfam_dicts = copy.deepcopy(lookup_uniprot_pfams[protein])  
except Exception as e:
    pfam_dicts = []  



def test_get_segments_to_model():  
    output = get_segments_to_model(protein, full_sequence, pfam_dicts, crystal_dicts)  
    assert output[0]["fasta_header"] == ">P16403_full_protein"
    assert output[1]["fasta_header"] == ">P16403_PF00538_Linker_histone_37_108"
    assert output[2]["fasta_header"] == ">P16403_PF00538_Linker_histone_1_213"
    assert len(output[0]["fasta_sequence"]) == len(full_sequence)
