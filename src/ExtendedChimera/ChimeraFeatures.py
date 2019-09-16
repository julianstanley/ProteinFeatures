# re for regular expression substitution
import re

# All of the features that we care about
features = ['AAindCodonDiv', 'AAindMolVol', 'AAindPolarity', 'AAindSS',
            "AAindElecCharge", "netCharge", "posCharge", "negCharge",
            "sumMetals", "CO", "CU", "FE", "K", "M", "MN", "MO", "NA",
            "NI", "ZN", "contacts", "RPKT", "Arg", "Lys", "Pro", "Thr",
            "surfMC", "surfM", "surfC", "hydro", "OHRxnConst", "SS", "helix",
            "sheet", "disord", "disordScore", "protBindSPPIDER",
            "protBindDISOPREDscore", "areaSAS", "reactivity",
            "circularVariance", "depth"]


def get_residue_features_metals(eResidue):
    residue_features = {
        "sumMetals_res": float('NaN'),
        "CO_res": eResidue.metals.count("CO"),
        "CU_res": eResidue.metals.count("CU"),
        "FE_res": eResidue.metals.count("FE"),
        "K_res": eResidue.metals.count("K"),
        "MG_res": eResidue.metals.count("MG"),
        "MN_res": eResidue.metals.count("MN"),
        "MO_res": eResidue.metals.count("MO"),
        "NA_res": eResidue.metals.count("NA"),
        "NI_res": eResidue.metals.count("NI"),
        "ZN_res": eResidue.metals.count("ZN")
    }


def get_residue_features(eResidue):
    ''' Get the features associated with this ExtendedResidue.
    '''
    residue_features = {
        "AA": eResidue.residue_1_letter,
        "Position": eResidue.number,
        "AAindCodonDiv_res": eResidue.aaind_codon_diversity,
        "AAindMolVol_res": eResidue.aaind_molecular_volume,
        "AAindPolarity_res": eResidue.aaind_polarity,
        "AAindSS_res": eResidue.aaind_secondary_structure,
        "AAindElecCharge_res": eResidue.aaind_charge,
        "netCharge_res": eResidue.charge,
        "posCharge_res": (1 if eResidue.charge == 1 else 0),
        "negCharge_res": (-1 if eResidue.charge == -1 else 0),
        'RPKT_res': (1 if eResidue.residue_1_letter in ["R", "K", "P", "T"]
                     else 0),
        "Arg_res": (1 if eResidue.residue_1_letter == "R" else 0),
        "Lys_res": (1 if eResidue.residue_1_letter == "K" else 0),
        "Pro_res": (1 if eResidue.residue_1_letter == "P" else 0),
        "Thr_res": (1 if eResidue.residue_1_letter == "T" else 0),
        "surfMC_res": float('NaN'),
        "surfM_res": float('NaN'),
        "surfC_res": float('NaN'),
        "hydro_res": eResidue.hydrophobicity,
        "OHRxnConst_res": eResidue.hydroxyl_constant,
        "SS_res": (1 if eResidue.isSheet or eResidue.isHelix else 0),
        "helix_res": (1 if eResidue.isHelix else 0),
        "sheet_res": (1 if eResidue.isSheet else 0),
        "disord_res": float('NaN'),
        "disordScore_res": float('NaN'),
        "protBindSPPIDER_res": float('NaN'),
        "protBindDISOPREDscore_res": float('NaN'),
        "areaSAS_res": eResidue.area_sas,
        "reactivity_res": eResidue.reactivity,
        "circularVariance_res": float('NaN'),
        "depth_res": eResidue.depth,
        "contacts_res": 1
    }

    # Make sure that we get all of the features
    for feature_name in features:
        res_feature_name = "{}_res".format(feature_name)
        if res_feature_name not in residue_features:
            residue_features[res_feature_name] = float('NaN')

    return residue_features


def get_residues_features_sums(eResidues):
    ''' Get the features associated with this list of ExtendedResidue
    '''
    residues_features = dict.fromkeys(features, 0)

    for residue in eResidues:
        # Update all sums with the features from this residue
        single_residue_features = get_residue_features(residue)
        for feature in residues_features:
            single_feature = "{}_res".format(feature)
            residues_features[feature] +=\
                single_residue_features[single_feature]

    return residues_features


def compute_global_attributes_residues(eResidues):
    """ Computes various residue-level attributes on the given residues

    Arguments:
        eResidues ([List of Extended Residues) A list of the residues
        on which to compute characteristics
    Requires:
        None
    Returns:
        0. sum_global_secondary_structure (int) How many residues have
        secondary structure?
        1. sum_global_helicies (int) How many residues are part of a helix?
        2. sum_global_sheets (int) How many residues are part of a sheet?
        3. sum_global_SAS (int)  How much solvent accessible surface area
        do all of the residues have in aggregate?
    Effect:
        None
    """
    return get_residues_features_sums(eResidues)


def compute_bubble_attributes_residues(base_atom, compared_residues, radius):
    base_atom.set_residue_contacts(compared_residues, radius)
    base_atom.set_metal_contacts(radius)
    return get_residues_features_sums(base_atom.residue_contacts)


def get_depth(all_atoms=False):
    """ Gets the depth of all RPKT to the surface

    Arguments:
    Requires:
    Returns:
    Effect:
    Note: This contains some repeated code, but it is probably necessary,
    as this collect distance data via runCommand instead of directly from
    atom objects. This is because I don't know how to measure the distance
    between an atom and the surface without runCommand.
    """
    # Read and clear the reply log before measuring distance
    # save_and_clear_reply_log(logfile_name)
    from ChimeraFeatures import *
    from ChimeraUtils import *
    from ExtendedAtom import *
    from Utils import *
    from chimera import selection

    if all_atoms:
        rc("select #0")
        rc("~select #0:@")
        rc(("measure distance :@ca"
            " selection multiple true show true"))
        r, distances = read_reply_log()

        rc("select protein")
        all_atoms = [ExtendedAtom(atom)
                     for atom in selection.currentAtoms()]

    else:
        rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys")
        all_atoms = [ExtendedAtom(atom)
                     for atom in selection.currentAtoms()]
        rc("select #0")
        rc("~select #0:@")
        rc(("measure distance :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys"
            " selection multiple true show true"))
        r, distances = read_reply_log()

        rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys")
        all_atoms = [ExtendedAtom(atom)
                     for atom in selection.currentAtoms()]

    # Map between atom number and atom residue
    atom_to_aa = {}
    for atom in all_atoms:
        atom_to_aa[atom.number] = atom.residue_1_letter

    # Initalize the return depth dictionary
    residue_to_depth = {}

    # Loop through each line that Chimera outputs
    for line in distances.splitlines():
        line = line.rstrip()
        if 'minimum distance from' in line:
            # Delete 'minmum distance from'
            line = line.replace('minimum distance from ', '')

            # Seperate the residue and distance by a tab
            line = re.sub(r" to #0:\? = ", "\t", line)

            # Get the atom and depth from the line
            atom, depth = line.split("\t")

            # Get just the residue number from the residue identifier
            atom_name = re.sub(r"#.+:", "", atom)
            atom_number = re.sub(r"@.+", "", atom_name)
            print(atom_to_aa)
            # Get the residue name to include the 1-letter aa identifier
            atom_with_aa = re.sub(
                r'@', (" " + atom_to_aa[atom_number] + " "), atom)

            residue_to_depth["{},{}".format(
                atom_number.split(".")[0], atom_to_aa[atom_number])] = depth

    print(residue_to_depth)

    return(residue_to_depth)


def get_depth_minimum(all_atoms=False):
    """ Gets the depth of all RPKT to the surface

    Arguments:
    Requires:
    Returns:
    Effect:
    Note: This contains some repeated code, but it is probably necessary,
    as this collect distance data via runCommand instead of directly from
    atom objects. This is because I don't know how to measure the distance
    between an atom and the surface without runCommand.
    """
    # Read and clear the reply log before measuring distance
    # save_and_clear_reply_log(logfile_name)
    from ChimeraFeatures import *
    from ChimeraUtils import *
    from ExtendedAtom import *
    from Utils import *
    from chimera import selection

    if all_atoms:
        rc("select #0")
        rc("~select #0:@")
        rc(("measure distance protein"
            " selection multiple true show true"))
        r, distances = read_reply_log()

        rc("select protein")
        all_atoms = [ExtendedAtom(atom)
                     for atom in selection.currentAtoms()]

    else:
        rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys")
        all_atoms = [ExtendedAtom(atom)
                     for atom in selection.currentAtoms()]
        rc("select #0")
        rc("~select #0:@")
        rc(("measure distance :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys"
            " selection multiple true show true"))
        r, distances = read_reply_log()

        rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys")
        all_atoms = [ExtendedAtom(atom)
                     for atom in selection.currentAtoms()]

    # Map between atom number and atom residue
    atom_to_aa = {}
    for atom in all_atoms:
        atom_to_aa[atom.number] = atom.residue_1_letter

    # Initalize the return depth dictionary
    residue_to_depth = {}

    # Loop through each line that Chimera outputs
    for line in distances.splitlines():
        line = line.rstrip()
        if 'minimum distance from' in line:
            # Delete 'minmum distance from'
            line = line.replace('minimum distance from ', '')

            # Seperate the residue and distance by a tab
            line = re.sub(r" to #0:\? = ", "\t", line)

            # Get the atom and depth from the line
            atom, depth = line.split("\t")

            residue_name, atom_type = atom.split("@")
            residue_name = re.sub(r"#.+:", "", atom)

            # Get just the residue number from the residue identifier
            atom_name = re.sub(r"#.+:", "", atom)
            atom_number = re.sub(r"@.+", "", atom_name)
            print(atom_to_aa)
            # Get the residue name to include the 1-letter aa identifier
            atom_with_aa = re.sub(
                r'@', (" " + atom_to_aa[atom_number] + " "), atom)

            residue_to_depth["{},{}".format(
                atom_number.split(".")[0], atom_to_aa[atom_number])] = depth

    print(residue_to_depth)

    return(residue_to_depth)
