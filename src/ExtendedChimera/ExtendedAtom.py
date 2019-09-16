from Utils import get_AAindex_df, charge_index_dict, aa_3_to_1, \
    reactivity_index_dict, OHrxnConst_index_dict, is_nonstandard_residue
import re
import chimera

# Import the AAindex dictionary from Utils
AAindex_df = get_AAindex_df()


class ExtendedResidue():
    ''' An ExtendedResidue (eResidue) has all of the properties
    of its first atom.
    '''

    def __init__(self, residue, depths=None, all_metals=[]):
        # Identifier attributes
        self.atoms = residue.atoms
        self.residue = residue
        self.residue_1_letter = aa_3_to_1[residue.type]
        self.name = str(residue)
        self.number = re.sub(r"#.+ .+ ", "", str(residue))  # .split(".")[0]
        print("Number here")
        print(self.number)

        # AAind attributes
        self.aaind_reactivity = reactivity_index_dict[self.residue_1_letter]
        self.aaind_charge = AAindex_df.loc[self.residue_1_letter,
                                           "AAindElecCharge"]
        self.aaind_codon_diversity = AAindex_df.loc[self.residue_1_letter,
                                                    "AAindCodonDiv"]
        self.aaind_molecular_volume = AAindex_df.loc[self.residue_1_letter,
                                                     "AAindMolVol"]
        self.aaind_secondary_structure = AAindex_df.loc[self.residue_1_letter,
                                                        "AAindSS"]
        self.aaind_polarity = AAindex_df.loc[self.residue_1_letter,
                                             "AAindPolarity"]

        # Residue structural attributes
        self.isSheet = residue.isSheet
        self.isHelix = residue.isHelix
        self.area_sas = residue.areaSAS

        # Misc. Attributes
        self.hydroxyl_constant = OHrxnConst_index_dict[self.residue_1_letter]
        self.charge = charge_index_dict[self.residue_1_letter]
        self.hydrophobicity = residue.kdHydrophobicity
        self.reactivity = reactivity_index_dict[self.residue_1_letter]

        # Depth
        if depths is None:
            self.depth = float("NaN")
            print("Depths is None for {}".format(self.name))
        else:
            if "{},{}".format(self.number.split(".")[0],
                              self.residue_1_letter) in depths:
                self.depth = depths["{},{}".format(self.number.split(".")[0],
                                                   self.residue_1_letter)]
            else:
                self.depth = float("NaN")
                print("Depth not found for {}".format(self.name))
                print("Number: {}".format(self.number))
                print("Number split: {}".format(self.number.split(".")))
                print("Tried key: {}").format(
                    "{},{}".format(self.number.split(".")[0],
                                   self.residue_1_letter))
                print("Depth dictionary: ")
                print(depths)

        # Direct metal contacts
        self.all_metals = all_metals
        direct_contacts = []
        for metal in all_metals:
            for atom in self.atoms:
                distance = chimera.distance(atom.xformCoord(),
                                            metal.metal_atom.xformCoord())
                if distance < metal.radius:
                    direct_contacts.append(metal)
                    break

        self.direct_metal_contacts = direct_contacts


class ExtendedAtom(ExtendedResidue, object):
    '''
    ExtendedAtom (eAtom)
    Define an extended extended atom class that holds a chimera atom and
    features of that atom. Holds an item as a property, instead of properly
    extending the atom class (since extension of the atom class is not allowed)
    '''

    def __init__(self, atom, depths=None, all_metals=[]):
        # This atom should have all of the attributes of its parent residue
        super(ExtendedAtom, self).__init__(atom.residue, depths, all_metals=[])

        self.residue = ExtendedResidue(atom.residue)

        # Identifier attributes
        self.atom = atom
        self.set_name()

        # The areaSAS of this atom
        self.set_area_sas(atom)

        # Set direct metal contacts
        self.set_metal_contacts(0)
        self.direct_metal_contacts = self.metal_contacts

    # Informal property setters
    def set_area_sas(self, atom):
        ''' Note: if surface is not computed before running this method,
        the area_sas will be 0 instead of None.
        '''
        try:
            self.area_sas = 0 if atom.areaSAS is None else atom.areaSAS
        except Exception:
            self.area_sas = 0

    def set_name(self):
        ''' Each ExtendedAtom has a name, formatted the way we like it
        '''
        # Get the residue associated with this atom
        residue = self.residue

        # Get the name of this atom, its number, and alternative locations
        atom_number = str(residue.residue)
        atom_small_name = self.atom.name
        alternative_loc = self.atom.altLoc
        alternative_loc = "" if alternative_loc == "" else "." \
            + alternative_loc

        # Initalize the name for the base atom
        atom_name = atom_number + " " + self.residue_1_letter + " "\
            + str(atom_small_name) + alternative_loc

        # Replace space + any three letters + space with ':'
        atom_name = re.sub(r"\ [A-Z][A-Z][A-Z]\ ", ":", atom_name)

        self.name = atom_name

    def set_atom_contacts(self, atoms, radius):
        ''' From the atoms given, returns a list of the atoms that are within
        a "radius" distance from this Extended Atom
        '''
        contacts = []
        for atom in atoms:
            if self.distance(atom) <= radius and self.distance(atom) >= 0:
                contacts.append(atom)
                print("here atom")
                print("The distance between {} and {} is {}").format(
                    self.name,
                    atom.name,
                    self.distance(atom))
        self.atom_contacts = contacts

    def set_residue_contacts(self, residues, radius):
        ''' For the residues given, returns a list of residues that are
        within a "radius" distance from this Extended Atom.
        '''
        contacts = []
        for residue in residues:
            # residue_eAtoms = [ExtendedAtom(atom) for atom in residue.atoms]
            # for residue_eAtom in residue_eAtoms:
            for res_atom in [ExtendedAtom(atom) for atom in residue.atoms]:
                # distance = self.distance(residue_eAtom)
                distance = self.distance(res_atom)
                if(distance <= radius and distance > 0 and not
                        is_nonstandard_residue(residue)):
                    contacts.append(residue)
                    # Break from this loop, since we only need
                    # to add the residue once
                    break

        self.residue_contacts = contacts

    def set_metal_contacts(self, radius):
        ''' For the metals given, sets the contacts of this residue that
        are within a "radius" distance from that metal
        '''
        metals = self.all_metals
        contacts = []
        for metal in metals:
            print('metal')
            print(metal)
            distance = self.distance(metal.metal_atom) - metal.radius
            if(distance <= radius):
                contacts.append(metal.metal_type)
        self.metal_contacts = contacts

    # Proper methods
    def distance(self, other_atom):
        ''' Gets the distance between this and another extendedAtom
        '''
        return(chimera.distance(self.atom.xformCoord(),
                                other_atom.atom.xformCoord()))

    def distances(self, list_other_eAtoms):
        '''
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
        '''
        distances = {}
        for eAtom in list_other_eAtoms:
            distances[eAtom.name] = self.distance(eAtom)
        return(distances)
