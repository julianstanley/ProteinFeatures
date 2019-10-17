from Utils import (
    get_AAindex_df,
    charge_index_dict,
    aa_3_to_1,
    reactivity_index_dict,
    OHrxnConst_index_dict,
    is_nonstandard_residue,
    ionicRadiusDict,
)
from ChimeraUtils import set_metal_contacts
import re
import chimera
import math

# Import the AAindex dictionary from Utils
AAindex_df = get_AAindex_df()


class ExtendedResidue:
    """ An ExtendedResidue (eResidue) has all of the properties
    of its first atom.
    """

    def __init__(
        self,
        residue,
        depths=None,
        all_metals=[],
        disopred_disorder=[],
        disopred_binding=[],
        sppider_binding=[],
    ):
        # Identifier attributes
        self.atoms = residue.atoms
        self.residue = residue
        try:
            self.residue_1_letter = aa_3_to_1[residue.type]
        # Better exception formatting
        except Exception as e:
            raise Exception(
                "{type} not found in residue {name}, {e}".format(
                    type=residue.type, name=str(residue)
                ),
                e=e,
            )
        self.name = str(residue)
        self.number = re.sub(r"#.+ .+ ", "", str(residue))  # .split(".")[0]

        # AAind attributes
        self.aaind_reactivity = reactivity_index_dict[self.residue_1_letter]
        self.aaind_charge = AAindex_df.loc[self.residue_1_letter, "AAindElecCharge"]
        self.aaind_codon_diversity = AAindex_df.loc[
            self.residue_1_letter, "AAindCodonDiv"
        ]
        self.aaind_molecular_volume = AAindex_df.loc[
            self.residue_1_letter, "AAindMolVol"
        ]
        self.aaind_secondary_structure = AAindex_df.loc[
            self.residue_1_letter, "AAindSS"
        ]
        self.aaind_polarity = AAindex_df.loc[self.residue_1_letter, "AAindPolarity"]

        # Residue structural attributes
        self.isSheet = residue.isSheet
        self.isHelix = residue.isHelix
        self.area_sas = residue.areaSAS

        # Misc. Attributes
        self.hydroxyl_constant = OHrxnConst_index_dict[self.residue_1_letter]
        self.charge = charge_index_dict[self.residue_1_letter]
        self.hydrophobicity = residue.kdHydrophobicity
        self.reactivity = reactivity_index_dict[self.residue_1_letter]

        # Disorder
        for score in disopred_disorder:
            if score["pos"] == self.number.split(".")[0]:
                self.disorder_score = score["Disorder Score"]
                self.disorder_call = score["Disorder Call"]
                if score["aa"] != self.residue_1_letter:
                    raise Exception(
                        "Mismatch: {}, {}, {}, {}, {}".format(
                            score,
                            self.number,
                            self.residue_1_letter,
                            residue.type,
                            aa_3_to_1,
                        )
                    )
                break

        # Disopred binding
        for d_binding in disopred_binding:
            if d_binding["pos"] == self.number.split(".")[0]:
                self.disopred_binding_score = d_binding["Binding Score"]
                self.disopred_binding_call = d_binding["Binding Call"]
                if score["aa"] != self.residue_1_letter:
                    raise Exception("Mismatch: {}, {}".format(d_binding, self.number))
                break

        # SPPIDER binding
        for s_binding in sppider_binding:
            if s_binding["pos"] == self.number.split(".")[0]:
                self.sppider_binding_call = s_binding["Binding Call"]
                if score["aa"] != self.residue_1_letter:
                    raise Exception("Mismatch: {}, {}".format(s_binding, self.number))
                break

        # Depth
        self.depths = depths
        if depths is None:
            self.depth = float("NaN")
        else:
            if (
                "{},{}".format(self.number.split(".")[0], self.residue_1_letter)
                in depths
            ):
                self.depth = depths[
                    "{},{}".format(self.number.split(".")[0], self.residue_1_letter)
                ]
            else:
                self.depth = float("NaN")

        # All metals in this structure
        self.all_metals = all_metals

        # Set inital contacts to be at a radius of 0
        self.metal_contacts = set_metal_contacts(self.atoms, all_metals, 0)

    def set_metal_contacts(self, radius):
        """ Set contacts between this residue and metals at a given radius
        """
        contacts = []

        for metal in self.all_metals:
            for atom in self.atoms:
                # Set the coordinates of this metal in a way that Point()
                # can accept
                metal_coords = [float(x) for x in metal.location.split(",")]
                # Get the distance between this atom and this metal
                distance = chimera.distance(
                    atom.xformCoord(),
                    chimera.Point(metal_coords[0], metal_coords[1], metal_coords[2]),
                )
                # If we're within the distance expected, add this metal
                # to the list and move to the next
                if distance - ionicRadiusDict[metal.type] <= radius:
                    contacts.append(metal)
                    break

        self.metal_contacts = contacts


class ExtendedAtom(ExtendedResidue, object):
    """
    ExtendedAtom (eAtom)
    Define an extended extended atom class that holds a chimera atom and
    features of that atom. Holds an item as a property, instead of properly
    extending the atom class (since extension of the atom class is not allowed)
    """

    def __init__(
        self,
        atom,
        depths=None,
        all_metals=[],
        disopred_disorder=[],
        disopred_binding=[],
        sppider_binding=[],
    ):
        # This atom should have all of the attributes of its parent residue
        super(ExtendedAtom, self).__init__(
            atom.residue,
            depths,
            all_metals,
            disopred_disorder,
            disopred_binding,
            sppider_binding,
        )

        # Identifier attributes
        self.atom = atom
        self.set_name()

        # The areaSAS of this atom
        self.set_area_sas(atom)

        # Set metal contacts at Inf to initalize
        self.set_metal_contacts(float("Inf"))
        self.global_metal_contacts = self.metal_contacts

    # Informal property setters
    def set_area_sas(self, atom):
        """ Note: if surface is not computed before running this method,
        the area_sas will be 0 instead of None.
        """
        try:
            self.area_sas = 0 if atom.areaSAS is None else atom.areaSAS
        except Exception:
            self.area_sas = 0

    def set_name(self):
        """ Each ExtendedAtom has a name, formatted the way we like it
        """

        # Get the name of this atom, its number, and alternative locations
        atom_number = str(self.residue)
        atom_small_name = self.atom.name
        alternative_loc = self.atom.altLoc
        alternative_loc = "" if alternative_loc == "" else "." + alternative_loc

        # Initalize the name for the base atom
        atom_name = (
            atom_number
            + " "
            + self.residue_1_letter
            + " "
            + str(atom_small_name)
            + alternative_loc
        )

        # Replace space + any three letters + space with ':'
        atom_name = re.sub(r"\ [A-Z][A-Z][A-Z]\ ", ":", atom_name)

        self.name = atom_name

    def set_atom_contacts(self, atoms, radius):
        """ From the atoms given, returns a list of the atoms that are within
        a "radius" distance from this Extended Atom. Does not include self.

        Takes in extended atoms
        """
        contacts = []
        for atom in atoms:
            if (
                self.distance(atom) <= radius
                and self.distance(atom) >= 0
                and atom.name != self.name
            ):
                contacts.append(atom)
        self.atom_contacts = contacts

    def set_residue_contacts(self, residues, radius):
        """ For the residues given, returns a list of residues that are
        within a "radius" distance from this Extended Atom.
        """
        contacts = []
        for residue in residues:
            # residue_eAtoms = [ExtendedAtom(atom) for atom in residue.atoms]
            # for residue_eAtom in residue_eAtoms:
            for res_atom in [ExtendedAtom(atom) for atom in residue.atoms]:
                # distance = self.distance(residue_eAtom)
                distance = self.distance(res_atom)
                if (
                    distance <= radius
                    and distance > 0
                    and not is_nonstandard_residue(residue)
                ):
                    contacts.append(residue)
                    # Break from this loop, since we only need
                    # to add the residue once
                    break

        self.residue_contacts = contacts

    def set_metal_contacts(self, radius):
        """ For the metals given, sets the contacts of this atom that
        are within a "radius" distance from that metal
        """
        metals = self.all_metals
        contacts = []
        for metal in metals:
            metal_coords = [float(x) for x in metal.location.split(",")]
            distance = chimera.distance(
                self.atom.xformCoord(),
                chimera.Point(metal_coords[0], metal_coords[1], metal_coords[2]),
            )
            if distance - ionicRadiusDict[metal.type] <= radius:
                contacts.append(metal)
        self.metal_contacts = contacts

    # Proper methods
    def distance(self, other_atom):
        """ Gets the distance between this and another extendedAtom
        """
        return chimera.distance(self.atom.xformCoord(), other_atom.atom.xformCoord())

    def distance_regular(self, other_atom):
        """ Gets the distance between this and a regular atom
        """
        return chimera.distance(self.atom.xformCoord(), other_atom.xformCoord())

    def distances(self, list_other_eAtoms):
        """
        Computes the iteratomic distances between this Extended Atom
        and the given list of Extended Atoms

        Arguments:
            eAtoms (List of ExtendedAtom): A list of all atoms to compare
            with the base atom.
        Requires:
            None
        Returns:
            (dict {eAtom --> float}): A mapping between each eAtom and the
            distance between that eAtom and the base eAtom.
        Effect:
            None
        """
        distances = {}
        for eAtom in list_other_eAtoms:
            distances[eAtom.name] = self.distance(eAtom)
        return distances

    def get_circular_variance(self):
        """ Gets the circular variance between this atom and all of it's
        current atomic neighbors
        """

        def unitVec(self_xformCoord, neighbor_xformCoord, index):
            """ Local function for defining unit vectors
            """
            return (neighbor_xformCoord[index] - self_xformCoord[index]) / math.sqrt(
                (neighbor_xformCoord[0] - self_xformCoord[0]) ** 2
                + (neighbor_xformCoord[1] - self_xformCoord[1]) ** 2
                + (neighbor_xformCoord[2] - self_xformCoord[2]) ** 2
            )

        # Accumulator vector
        sumUnitVectors = [0, 0, 0]

        # Get the vector for all neighbors
        neighbors = self.atom_contacts
        for neighbor in neighbors:
            self_xformCoord = [
                float(coord) for coord in str(self.atom.xformCoord()).split(" ")
            ]

            neighbor_xformCoord = [
                float(coord) for coord in str(neighbor.atom.xformCoord()).split(" ")
            ]

            unitVecAB_x = unitVec(self_xformCoord, neighbor_xformCoord, 0)
            unitVecAB_y = unitVec(self_xformCoord, neighbor_xformCoord, 1)
            unitVecAB_z = unitVec(self_xformCoord, neighbor_xformCoord, 2)

            # Accumulate the vector with this neighbor
            sumUnitVectors = [
                sumUnitVectors[0] + unitVecAB_x,
                sumUnitVectors[1] + unitVecAB_y,
                sumUnitVectors[2] + unitVecAB_z,
            ]

        magSumUnitVectors = math.sqrt(
            (sumUnitVectors[0]) ** 2
            + (sumUnitVectors[1]) ** 2
            + (sumUnitVectors[2]) ** 2
        )

        # The circular variance is 1 - the unit vector magnitude, normalized to count
        # try-except is necessary for when bubble attributes are calculated at a radius
        # of zero.
        try:
            cv = 1 - (magSumUnitVectors / len(neighbors))
        except Exception as e:
            print(e)
            cv = 0

        return cv
