import pandas as pd

reactivity_index_dict = {
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 1.009211433,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 1.296184165,
    "Q": 0,
    "R": 0.900412937,
    "S": 0,
    "T": 0.830171969,
    "V": 0,
    "W": 0,
    "Y": 0,
    "X": 0,
}

charge_index_dict = {
    "A": 0,
    "C": 0,
    "D": -1,
    "E": -1,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 1,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 1,
    "S": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0,
    "X": 0,
}

OHrxnConst_index_dict = {
    "A": 77000000,
    "C": 34000000000,
    "D": 75000000,
    "E": 230000000,
    "F": 6500000000,
    "G": 17000000,
    "H": 13000000000,
    "I": 1800000000,
    "K": 340000000,
    "L": 1700000000,
    "M": 8300000000,
    "N": 49000000,
    "P": 480000000,
    "Q": 540000000,
    "R": 3500000000,
    "S": 320000000,
    "T": 510000000,
    "V": 760000000,
    "W": 13000000000,
    "Y": 13000000000,
}

kd_hydrophobicity = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

aa_3_to_1 = aaDict = dict(
    {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLN": "Q",
        "GLU": "E",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "PTR": "Y",
        "CSX": "C",
        "CSO": "C",
        "CSD": "C",
        "KCX": "K",
        "UNK": "X",
        "MLY": "K",
        "PVH": "H",
    }
)

ionicRadiusDict = dict(
    {
        "CA": 1.00,
        "CO": 0.65,
        "CU": 0.73,
        "FE": 0.61,
        "K": 1.38,
        "MG": 0.72,
        "MN": 0.83,
        "MO": 0.65,
        "NA": 1.02,
        "NI": 0.69,
        "ZN": 0.74,
    }
)


def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n // 2 - 1 : n // 2 + 1]) / 2.0, s[n // 2])[n % 2] if n else None


def is_nonstandard_residue(residue):
    return (
        "MSE" in str(residue)
        or "CSX" in str(residue)
        or "UNK" in str(residue)
        or "CSD" in str(residue)
        or "CSO" in str(residue)
        or "KCX" in str(residue)
    )


def make_standard(
    protein, constant_features, atomic_features, atom_names, check_atomic_features=True
):
    """
  Format:
  Keys: (Protein Name, Atom Name, Feature Name)
  Values: Value of this feature for this atom in this protein
  """
    features = {}

    # Check the atomic features to make sure all of their keys are in atom_name
    # This is resource-heavy, but a nice-to-have
    if check_atomic_features:
        for atomic_feature in atomic_features:
            for atom_name in atom_names:
                try:
                    atomic_feature["Value"][atom_name]
                except Exception:
                    raise Exception(
                        """
              The atomic feature {} doesn't have an
              entry for the atom {}.
              """.format(
                            atomic_feature["Name"], atom_name
                        )
                    )
                    print("Checked atomic features")

    # We need an entry for each atom name
    for atom_name in atom_names:
        # Add entries for the constant features
        for constant_feature in constant_features:
            feature_name = constant_feature["Name"]
            feature_value = constant_feature["Value"]
            features[(protein, atom_name, feature_name)] = feature_value

        # Add entries for the atomic features
        for atomic_feature in atomic_features:
            feature_name = atomic_feature["Name"]
            atom_to_feature = atomic_feature["Value"]
            feature_value = atom_to_feature[atom_name]
            features[(protein, atom_name, feature_name)] = feature_value

            return features


def get_AAindex_df():
    """
  Returns: (Pandas DataFrame Object) - Columns represent AAindex
  characteristics, rows represent amino acids
  """
    AAindex_mat = [
        [-0.591, -1.302, -0.733, 1.57, -0.146],
        [-1.343, 0.465, -0.862, -1.020, -0.255],
        [1.05, 0.302, -3.656, -0.259, -3.242],
        [1.357, -1.453, 1.477, 0.113, -0.837],
        [-1.006, -0.590, 1.891, -0.397, 0.412],
        [-0.384, 1.652, 1.33, 1.045, 2.064],
        [0.336, -0.417, -1.673, -1.474, -0.078],
        [-1.239, -0.547, 2.131, 0.393, 0.816],
        [1.831, -0.561, 0.533, -0.277, 1.648],
        [-1.019, -0.987, -1.505, 1.266, -0.912],
        [-0.663, -1.524, 2.219, -1.005, 1.212],
        [0.945, 0.828, 1.299, -0.169, 0.933],
        [0.189, 2.081, -1.628, 0.421, -1.392],
        [0.931, -0.179, -3.005, -0.503, -1.853],
        [1.538, -0.055, 1.502, 0.44, 2.897],
        [-0.228, 1.399, -4.760, 0.67, -2.647],
        [-0.032, 0.326, 2.213, 0.908, 1.313],
        [-1.337, -0.279, -0.544, 1.242, -1.262],
        [-0.595, 0.009, 0.672, -2.128, -0.184],
        [0.26, 0.83, 3.097, -0.838, 1.512],
    ]

    AAindex_df = pd.DataFrame(AAindex_mat)
    AAindex_df["AAs"] = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]
    AAindex_df.set_index("AAs", inplace=True)
    AAindex_df.columns = [
        "AAindPolarity",
        "AAindSS",
        "AAindMolVol",
        "AAindCodonDiv",
        "AAindElecCharge",
    ]
    return AAindex_df