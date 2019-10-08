# re for regular expression substitution
import re
from ChimeraUtils import read_reply_log
from chimera import selection
from chimera import runCommand as rc
from ExtendedAtom import ExtendedAtom
from Utils import median

# All of the features that we care about
features = ['AAindCodonDiv', 'AAindMolVol', 'AAindPolarity', 'AAindSS',
            "AAindElecCharge", "netCharge", "posCharge", "negCharge",
            "sumMetals", "CO", "CU", "FE", "K", "M", "MN", "MO", "NA",
            "NI", "ZN", "contacts", "RPKT", "Arg", "Lys", "Pro", "Thr",
            "surfMC", "surfM", "surfC", "hydro", "OHRxnConst", "SS", "helix",
            "sheet", "disord", "disordScore", "protBindSPPIDER",
            "protBindDISOPREDscore", "areaSAS", "reactivity",
            "circularVariance", "depth"]


def get_atom_features(eAtom):
    ''' Get the features associated with this ExtendedAtom
    '''
    atom_features = {
        "AA": eAtom.residue_1_letter,
        "Position": eAtom.number,
        "AAindCodonDiv": eAtom.aaind_codon_diversity,
        "AAindMolVol": eAtom.aaind_molecular_volume,
        "AAindPolarity": eAtom.aaind_polarity,
        "AAindSS": eAtom.aaind_secondary_structure,
        "AAindElecCharge": eAtom.aaind_charge,
        "netCharge": eAtom.charge,
        "posCharge": (1 if eAtom.charge == 1 else 0),
        "negCharge": (-1 if eAtom.charge == -1 else 0),
        # Keep sumMetals as NaN, since metals sum isn't strictly addative
        "sumMetals": len(list(set(
            [metal.site_number for metal in eAtom.metal_contacts]))),
        "CA": eAtom.metal_contacts.count("CA"),
        "CO": eAtom.metal_contacts.count("CO"),
        "CU": eAtom.metal_contacts.count("CU"),
        "FE": eAtom.metal_contacts.count("FE"),
        "K": eAtom.metal_contacts.count("K"),
        "MG": eAtom.metal_contacts.count("MG"),
        "MN": eAtom.metal_contacts.count("MN"),
        "MO": eAtom.metal_contacts.count("MO"),
        "NA": eAtom.metal_contacts.count("NA"),
        "NI": eAtom.metal_contacts.count("NI"),
        "ZN": eAtom.metal_contacts.count("ZN"),
        'RPKT': (1 if eAtom.residue_1_letter in ["R", "K", "P", "T"] else 0),
        "Arg": (1 if eAtom.residue_1_letter == "R" else 0),
        "Lys": (1 if eAtom.residue_1_letter == "K" else 0),
        "Pro": (1 if eAtom.residue_1_letter == "P" else 0),
        "Thr": (1 if eAtom.residue_1_letter == "T" else 0),
        "surfMC": float('NaN'),
        "surfM": float('NaN'),
        "surfC": float('NaN'),
        "hydro": eAtom.hydrophobicity,
        "OHRxnConst": eAtom.hydroxyl_constant,
        "SS": (1 if eAtom.isSheet or eAtom.isHelix else 0),
        "helix": (1 if eAtom.isHelix else 0),
        "sheet": (1 if eAtom.isSheet else 0),
        "disord": float('NaN'),
        "disordScore": float('NaN'),
        "protBindSPPIDER": float('NaN'),
        "protBindDISOPREDscore": float('NaN'),
        "areaSAS": eAtom.area_sas,
        "reactivity": eAtom.reactivity,
        "circularVariance": float('NaN'),
        "depth": eAtom.depth,
        "contacts": 1
    }

    # Make sure that we get all of the features
    for feature_name in features:
        if feature_name not in atom_features:
            atom_features[feature_name] = float('NaN')

    return atom_features


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
        # Keep sumMetals as NaN, since metals sum isn't strictly addative
        "sumMetals_res": len(list(set(
            [metal.site_number for metal in eResidue.metal_contacts]))),
        "CA_res": eResidue.metal_contacts.count("CA"),
        "CO_res": eResidue.metal_contacts.count("CO"),
        "CU_res": eResidue.metal_contacts.count("CU"),
        "FE_res": eResidue.metal_contacts.count("FE"),
        "K_res": eResidue.metal_contacts.count("K"),
        "MG_res": eResidue.metal_contacts.count("MG"),
        "MN_res": eResidue.metal_contacts.count("MN"),
        "MO_res": eResidue.metal_contacts.count("MO"),
        "NA_res": eResidue.metal_contacts.count("NA"),
        "NI_res": eResidue.metal_contacts.count("NI"),
        "ZN_res": eResidue.metal_contacts.count("ZN"),
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
    ''' Get the sum of all the features for the residues in the given
    list of ExtendedResidues.
    '''
    residues_features = dict.fromkeys(features, 0)

    for residue in eResidues:
        # Get the features for this single residue
        single_residue_features = get_residue_features(residue)
        # Add this residue's features to the accumulating features
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

    features = get_residues_features_sums(eResidues)

    # Overwrite metal binding features
    # Metal binding features at a global level will not vary residue-to-residue,
    # so just get all of the metals infinite distance from one of the residues
    eResidue = eResidues[0]
    eResidue.set_metal_contacts(float('Inf'))
    metal_types = [metal.type for metal in eResidue.metal_contacts]

    # Update the features dictionary with contacts from this residue
    features.update({
        "sumMetals": len(
            list(set([metal.site_number for metal in eResidue.metal_contacts]))),
        "CA": metal_types.count("CA"),
        "CO": metal_types.count("CO"),
        "CU": metal_types.count("CU"),
        "FE": metal_types.count("FE"),
        "K": metal_types.count("K"),
        "MG": metal_types.count("MG"),
        "MN": metal_types.count("MN"),
        "MO": metal_types.count("MO"),
        "NA": metal_types.count("NA"),
        "NI": metal_types.count("NI"),
        "ZN": metal_types.count("ZN")})

    return(features)


def compute_bubble_attributes_residues(base_atom, compared_residues, radius):
    base_atom.set_residue_contacts(compared_residues, radius)
    base_atom.set_metal_contacts(radius)
    bubble_features = get_residues_features_sums(base_atom.residue_contacts)

    # Add the base features to the bubble features, if we're at a radius of 0
    if radius < 1e-20:
        base_features = get_atom_features(base_atom)
        for key in bubble_features:
            # Don't include self in measurement of RPKT and contacts
            if key in base_features and key not in ["R", "K", "P", "T", "RPKT",
                                                    "contacts"]:
                bubble_features[key] = bubble_features[key] + base_features[key]

    # Overwrite metal binding features
    # Metal binding features at a bubble level does not depend on compared residues,
    # just the base atom
    base_atom.set_metal_contacts(radius)
    metal_types = [metal.type for metal in base_atom.metal_contacts]

    # Update the features dictionary with contacts from this residue
    bubble_features.update({
        "sumMetals": len(
            list(set([metal.site_number for metal in base_atom.metal_contacts]))),
        "CA": metal_types.count("CA"),
        "CO": metal_types.count("CO"),
        "CU": metal_types.count("CU"),
        "FE": metal_types.count("FE"),
        "K": metal_types.count("K"),
        "MG": metal_types.count("MG"),
        "MN": metal_types.count("MN"),
        "MO": metal_types.count("MO"),
        "NA": metal_types.count("NA"),
        "NI": metal_types.count("NI"),
        "ZN": metal_types.count("ZN")})

    return(bubble_features)


def get_depths(times=5):
    with open("troubleshoot.txt", 'a') as file:
        file.write("In get_depths")

    depths_lists = {}
    depths_return = {}
    for time in range(times):
        current_depths = get_depth()
        for atom_name, depth in current_depths.items():
            if atom_name in depths_lists:
                depths_lists[atom_name].append(depth)
            else:
                depths_lists[atom_name] = [depth]

    with open("troubleshoot.txt", 'a') as file:
        file.write("First loop")

    for atom_name, depths_list in depths_lists.items():
        depths_return[atom_name] = median(
            [float(depth) for depth in depths_list])

    with open("troubleshoot.txt", 'a') as file:
        file.write("Second loop")

    return(depths_return)


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

    with open("troubleshoot.txt", 'a') as file:
        file.write("Get_depth")

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
        with open("troubleshoot.txt", 'a') as file:
            file.write("Before rc")
        rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys")
        all_atoms = [ExtendedAtom(atom)
                     for atom in selection.currentAtoms()]
        rc("select #0")
        rc("~select #0:@")
        rc(("measure distance :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys"
            " selection multiple true show true"))
        r, distances = read_reply_log()

        rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb|:met|:cys")
        with open("troubleshoot.txt", 'a') as file:
            file.write("Before init atoms")
        all_atoms = [ExtendedAtom(atom)
                     for atom in selection.currentAtoms()]
        with open("troubleshoot.txt", 'a') as file:
            file.write("After init atoms")

    with open("troubleshoot.txt", 'a') as file:
        file.write("Get_depth after first block")

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
