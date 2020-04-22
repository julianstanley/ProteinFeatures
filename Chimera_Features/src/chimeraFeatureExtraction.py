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
    attempts_limit=5,
):
    """
    Arguments:
    metal_binding_sites: List of metal binding sites, each in dictionary
    keys: 'File Name', 'Site Number', 'Metal ID', 'Residues'.
    Global Features:
    """
    # Open the pdb file, generate a surface, then perform DockPrep
    # (Add hydrogens, etc).

    # Try, in case surface generation fails.
    attempts = 0
    while attempts < attempts_limit:
        try:
            rc("open " + pdb_file)
            generate_surface()
            break
        except Exception as e:
            attempts += 1
            if attempts >= attempts_limit:
                raise Exception(
                    "Surface computation failed {} times".format(attempts))

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
        coords = re.sub(
            r"centroid name, ID, center: centroid: c1 \( *", "", coords)
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

    # Get RPKT atoms and residues
    rc("~select")
    rc("select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb")
    print("getting all {} rpkt residues".format(
        len(selection.currentResidues())))

    rpkt_residues = [
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

    # Get atoms and residues near RPKT
    # Efficiency technique: include already-calculated RPKT atoms, delete those from selection
    rc("~select")
    rc("select protein")
    rc(
        "select :arg@cd z<12|:lys@ce z<12|:mly@ce z<12|:kcx@ce z<12|:pro@cd z<12|:thr@cb z<12"
    )

    residues_near_rpkt = list()
    for residue in selection.currentResidues():
        residues_near_rpkt.append(
            ExtendedResidue(
                residue,
                depths,
                metals,
                disopred_disorder_map,
                disopred_binding_map,
                sppider_binding_map,
            )
        )

    atoms_near_rpkt = list()
    for atom in selection.currentAtoms():
        atoms_near_rpkt.append(
            ExtendedAtom(
                atom,
                depths,
                metals,
                disopred_disorder_map,
                disopred_binding_map,
                sppider_binding_map,
            )
        )

    # Get all atoms and residues
    # Efficiency technique: include already-calculated atoms near RPKT, delete those from selection
    rc("~select")
    rc("select protein")
    rc(
        "~select :arg@cd z<12|:lys@ce z<12|:mly@ce z<12|:kcx@ce z<12|:pro@cd z<12|:thr@cb z<12"
    )
    rc("~select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb")
    print("getting all {} additional residues".format(
        len(selection.currentResidues())))
    all_residues = list(residues_near_rpkt)
    for residue in selection.currentResidues():
        all_residues.append(
            ExtendedResidue(
                residue,
                depths,
                metals,
                disopred_disorder_map,
                disopred_binding_map,
                sppider_binding_map,
            )
        )

    # No need to get all atoms within 12A.
    print("done")

    print("getting global attributes")
    global_attributes = compute_global_attributes_residues(all_residues)
    global_attributes.pop("circularVariance")

    save_and_clear_reply_log(logfile)

    print("getting bubble attributes")
    atom_radius_features = {}
    for radius in radii:
        radius_plus_bit = radius + 1e-20
        print("Analyzing at radius: {}\n".format(radius))
        # Get atoms and residues near RPKT atoms
        rc("~select")
        rc(
            (
                "select :arg@cd z<{radius}|:lys@ce z<{radius}|"
                ":mly@ce z<{radius}| :kcx@ce z<{radius}|:pro@cd z<{radius}|"
                ":thr@cb z<{radius}"
            ).format(radius=radius_plus_bit)
        )
        rc("~select :arg@cd|:lys@ce|:mly@ce|:kcx@ce|:pro@cd|:thr@cb")
        print(
            "getting all {} non-RPKT residues near rpkt".format(
                len(selection.currentResidues())
            )
        )

        # Get all residues in the selection (near RPKT, building on previous RPKT list)
        residues_near_rpkt = list(rpkt_residues)
        for residue in selection.currentResidues():
            residues_near_rpkt.append(
                ExtendedResidue(
                    residue,
                    depths,
                    metals,
                    disopred_disorder_map,
                    disopred_binding_map,
                    sppider_binding_map,
                )
            )

        atoms_near_rpkt = list(rpkt_atoms)
        for atom in selection.currentAtoms():
            atoms_near_rpkt.append(
                ExtendedAtom(
                    atom,
                    depths,
                    metals,
                    disopred_disorder_map,
                    disopred_binding_map,
                    sppider_binding_map,
                )
            )

        print("done")

        # Get bubble attributes
        for atom in rpkt_atoms:
            if atom.name in atom_radius_features:
                atom_radius_features[atom.name][
                    radius
                ] = compute_bubble_attributes_residues(
                    atom, all_residues, atoms_near_rpkt, radius
                )
            else:
                # This must be the first radius the loop has seen.
                # atom_radius_features[atom.name] should point to a dict
                # with radii for keys, bubble attributes for values
                atom_radius_features[atom.name] = {
                    radius: compute_bubble_attributes_residues(
                        atom, all_residues, atoms_near_rpkt, radius
                    )
                }

    print("getting atom-level attributes")

    atom_features = {}
    for atom in rpkt_atoms:
        eResidue = ExtendedResidue(atom.residue, depths, metals)
        extended_atoms = []
        for a in eResidue.atoms:
            extended_atoms.append(ExtendedAtom(a, depths, metals))

        atom_features[atom.name] = compute_bubble_attributes_residues(
            atom, residues_near_rpkt, extended_atoms, 0
        )

        # Format each feature with _res to match labels of original pipeline data
        for key in atom_features[atom.name]:
            if "res" not in key:
                atom_features[atom.name]["{}_res".format(key)] = atom_features[
                    atom.name
                ].pop(key)

        # Add global circular variance
        atom_features[atom.name].update(
            compute_global_circular_variance(atom, atoms_near_rpkt)
        )

        atom_features[atom.name].update(
            {"AA": eResidue.residue_1_letter, "Position": eResidue.number}
        )

    return (global_attributes, atom_features, atom_radius_features)
