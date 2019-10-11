from ExtendedChimera.ChimeraFeatures import (
    get_depths,
    compute_bubble_attributes_residues,
    compute_global_attributes_residues,
    compute_global_circular_variance,
)
from ExtendedChimera.MetalAtom import MetalAtom
from ExtendedChimera.ChimeraUtils import (
    generate_surface,
    add_hydrogens_prep,
    save_and_clear_reply_log,
    read_reply_log,
)
from ExtendedChimera.ExtendedAtom import ExtendedAtom, ExtendedResidue
from ExtendedChimera.Utils import ionicRadiusDict
from chimera import selection
from chimera import runCommand as rc
import re


def chimeraFeatureExtraction(
    pdb_file,
    radii=[5, 8, 10, 12],
    metal_binding_sites=[],
    disopred_disorder_map=[],
    disopred_binding_map=[],
    sppider_binding_map=[],
    logfile="log.txt",
):
    """
    Arguments:
    metal_binding_sites: List of metal binding sites, each in dictionary
    keys: 'File Name', 'Site Number', 'Metal ID', 'Residues'.
    Global Features:
    """
    # Open the pdb file, generate a surface, then perform DockPrep
    # (Add hydrogens, etc).
    rc("open " + pdb_file)
    generate_surface()
    rc("~select")
    add_hydrogens_prep()

    # Initalize metal binding features for this pdb file
    metals = []
    for metal_binding_site in metal_binding_sites:
        metal_type = metal_binding_site["Metal ID"]
        site_number = metal_binding_site["Site Number"]
        # Convert residues from:
        # C208,C211,etc ; to:
        # #:208@|#:211@
        residues = metal_binding_site["Residues"]
        residues_selection = re.sub(r",", "@|#:", residues)
        residues_selection = re.sub(r"[a-zA-Z]", "", residues_selection)
        residues_selection = "#:{}@".format(residues_selection)

        # Select all residues associated with this metal binding site
        rc("select {}".format(residues_selection))

        # Deselect backbone (only side chains should be selected)
        rc("~select @n,ca,c,o")

        # Clear the reply log, since we'll need it in a minute
        save_and_clear_reply_log(logfile)

        # Define a centroid around the currently selected residues
        rc(
            (
                "define centroid massWeighting false radius {} raiseTool false "
                "number 1 selection"
            ).format(ionicRadiusDict[metal_type])
        )

        # Get the x,y,z coordinates of the metal from the reply log
        r, coords = read_reply_log()
        coords = re.sub(r"centroid name, ID, center: centroid: c1 \( *", "", coords)
        coords = re.sub(r"\)", "", coords)
        coords = re.sub(r",", "", coords)
        coords = re.sub(r" +", " ", coords)
        coords = re.sub(r"\n", "", coords)
        coords = re.sub(r"\r", "", coords)
        coords = ",".join(coords.split())
        save_and_clear_reply_log(logfile)

        metals.append(MetalAtom(metal_type, residues, coords, site_number))

    # Depths
    depths = get_depths()

    # New
    save_and_clear_reply_log(logfile)

    # Get RPKT atoms and residues
    rc("~select")
    rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb")
    # rpkt_residues = [ExtendedResidue(residue)
    #                 for residue in selection.currentResidues()]
    print("getting all {} rpkt residues".format(len(selection.currentResidues())))
    # New
    save_and_clear_reply_log(logfile)
    rpkt_atoms = [
        ExtendedAtom(
            atom,
            depths,
            metals,
            disopred_disorder_map,
            disopred_binding_map,
            sppider_binding_map,
        )
        for atom in selection.currentAtoms()
    ]
    print("done")
    save_and_clear_reply_log(logfile)

    # Get all atoms and residues
    rc("~select")
    rc("select protein")
    print("getting all {} residues".format(len(selection.currentResidues())))
    all_residues = [
        ExtendedResidue(
            residue,
            depths,
            metals,
            disopred_disorder_map,
            disopred_binding_map,
            sppider_binding_map,
        )
        for residue in selection.currentResidues()
    ]

    print("done")
    atoms_near_rpkt = [
        ExtendedAtom(
            atom,
            depths,
            metals,
            disopred_disorder_map,
            disopred_binding_map,
            sppider_binding_map,
        )
        for atom in selection.currentAtoms()
    ]

    print("getting global attributes")
    global_attributes = compute_global_attributes_residues(all_residues)
    global_attributes.pop("circularVariance")

    save_and_clear_reply_log(logfile)

    print("getting bubble attributes")
    atom_radius_features = {}
    for radius in radii:
        radius_plus_bit = radius + 1e-20
        print("Analyzing at radius: {}".format(radius))
        # Get atoms and residues near RPKT atoms
        rc("~select")
        rc(
            (
                "select :arg@cd z<{radius}|:lys@ce z<{radius}|"
                ":mly@ce z<{radius}| :kcx@ce z<{radius}|:pro@cd z<{radius}|"
                ":thr@cb z<{radius}"
            ).format(radius=radius_plus_bit)
        )
        print(
            "getting all {} residues near rpkt".format(len(selection.currentResidues()))
        )

        # Added try-except due to error on 09/18
        residues_near_rpkt = []
        for residue in selection.currentResidues():
            try:
                eResidue = ExtendedResidue(
                    residue,
                    depths,
                    metals,
                    disopred_disorder_map,
                    disopred_binding_map,
                    sppider_binding_map,
                )
                residues_near_rpkt.append(eResidue)
            except Exception:
                continue

        print("done")

        # Get bubble attributes
        for atom in rpkt_atoms:
            if atom.name in atom_radius_features:
                atom_radius_features[atom.name][
                    radius
                ] = compute_bubble_attributes_residues(
                    atom, residues_near_rpkt, atoms_near_rpkt, radius
                )
            else:
                # This must be the first radius the loop has seen.
                # atom_radius_features[atom.name] should point to a dict
                # with radii for keys, bubble attributes for values
                atom_radius_features[atom.name] = {
                    radius: compute_bubble_attributes_residues(
                        atom, residues_near_rpkt, atoms_near_rpkt, radius
                    )
                }

    print("getting atom-level attributes")
    all_atoms = [
        ExtendedAtom(
            atom,
            depths,
            metals,
            disopred_disorder_map,
            disopred_binding_map,
            sppider_binding_map,
        )
        for atom in selection.currentAtoms()
    ]
    atom_features = {}
    for atom in rpkt_atoms:
        eResidue = ExtendedResidue(atom.residue, depths, metals)
        atom_features[atom.name] = compute_bubble_attributes_residues(
            atom,
            all_residues,
            [ExtendedAtom(a, depths, metals) for a in eResidue.atoms],
            0,
        )

        # Format each feature with _res to match labels of original pipeline data
        for key in atom_features[atom.name]:
            if "res" not in key:
                atom_features[atom.name]["{}_res".format(key)] = atom_features[
                    atom.name
                ].pop(key)

        # Add global circular variance
        atom_features[atom.name].update(
            compute_global_circular_variance(atom, all_atoms)
        )

        atom_features[atom.name].update(
            {"AA": eResidue.residue_1_letter, "Position": eResidue.number}
        )

    return (global_attributes, atom_features, atom_radius_features)
