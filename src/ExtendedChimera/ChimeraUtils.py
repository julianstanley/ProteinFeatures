from chimera import runCommand as rc
from chimera import dialogs
from DockPrep import prep
import chimera
import AddH
from DockPrep.prefs import defaults, INCOMPLETE_SC
from Utils import ionicRadiusDict


def read_reply_log():
    """ Reads the Chimera reply log

    Arguments:
        None
    Requires:
        There must be a Chimera reply log open. This means that
        Chimera needs to be running in GUI mode, e.g. Chimera needs to be
        opened with the command `chimera --start "Reply Log" ...`
    Returns:
        r (Chimera dialog object): The reply log object
        r_text (str): The current reply log text
    Effect:
        None
    """
    r = dialogs.find('reply')
    r_text = r.text.get('1.0', 'end')
    return (r, r_text)


def save_and_clear_reply_log(logfile_name):
    """ Reads the reply log, saves it in a file with the given name,
    and clears the reply log's contents

    Arguments:
        logfile_name (str): The name of the logfile to which the reply log
        should be saved
    Requires:
        There must be a Chimera reply log open. This means that
        Chimera needs to be running in GUI mode, e.g. Chimera needs to be
        opened with the command `chimera --start "Reply Log" ...`
    Returns:
        r (Chimera dialog object): The reply log object
    Effect:
        Clears any text from the reply log object
    """
    r, r_text = read_reply_log()
    with open(logfile_name, 'a') as logfile:
        logfile.write(r_text)
    r.Clear()
    return(r)


def clean_protein():
    """ Cleans the currently-loaded protein in Chimera

    Arguments:
        None
    Requires:
        A protein must be opened in Chimera (
        e.g. 'runCommand("open " + file_name))' must have already run
    Returns:
        Void
    Effect:
        Deletes the following items from the loaded protein structure:
        1. Non-peptide atoms
        2. ligands
        3. alternative residue locations
    """
    # Clean out all non-peptide atoms and ligands
    rc("select protein")
    rc("select invert")
    rc("delete selected")
    rc("select ligand")
    rc("delete selected")

    # Delete all alternative locations for residues
    rc("select #0:@.B")
    rc("delete selected")


def generate_surface():
    """ Generates a surface onto the currently-loaded protein in Chimera.
    If a surface cannot be generated, return an exception.
    Arguments:
        None
    Requires:
        A protein must be opened in Chimera (
        e.g. 'runCommand("open " + file_name))' must have already run
        It is also recommended that the loaded protein has already been
        cleaned with clean_protein()
    Returns:
        Void
    Effect:
        Attempts to split the current protein and then generates a
        surface for each part of the protein
    """
    # Select and split the protein
    rc("select protein")
    rc("split")

    # Generate the surface
    # Note: allComponents false excludes bubbles from surface computation
    rc("surface allComponents false")

    # Check surface for computation failure
    r, r_text = read_reply_log()
    if 'connected surface components' not in r_text:
        rc("close all")
        raise Exception("Surface computation failed")


def add_hydrogens_prep():
    """ Performs DockPrep on the currently-loaded protein in Chimera
    Arguments:
        None
    Requires:
        A protein must be opened in Chimera (
        e.g. 'runCommand("open " + file_name))' must have already run
        It is also recommended that the loaded protein has already been
        cleaned with clean_protein()
        If surface computation is run after this function, the surface
        computation may fail due to the added hydrogens
    Returns:
        Void
    Effect:
        Modifies the structure in the following ways:
        1. Adds hydrogens
        2. Mutates non-standard AAs to standard AAs
        3. Modifies incomplete side chains
        4. Deletes solven molecules
    """
    # Add hydrogen atoms/charges, mutate MSE and other non-standard
    # AAs to standard, and perform other basic DockPrep operations.
    models = chimera.openModels.list(modelTypes=[chimera.Molecule])
    # try:
    prep(models, addHFunc=AddH.hbondAddHydrogens, hisScheme=None,
         mutateMSE=True,
         mutate5BU=True, mutateUMS=True, mutateCSL=True, delSolvent=True,
         delIons=False, delLigands=False, delAltLocs=True,
         incompleteSideChains="rotamers", nogui=True,
         rotamerLib=defaults[INCOMPLETE_SC], rotamerPreserve=True,
         memorize=False, memorizeName=None)
    # except Exception as e:
    #   print(e)


def set_metal_contacts(atoms, all_metals, radius):
    ''' Set contacts between the given atoms and metals at a given radius
    '''
    contacts = []

    for metal in all_metals:
        for atom in atoms:
            # Set the coordinates of this metal in a way that Point()
            # can accept
            metal_coords = [float(x)
                            for x in metal.location.split(",")]
            # Get the distance between this atom and this metal
            distance = chimera.distance(atom.xformCoord(),
                                        chimera.Point(metal_coords[0],
                                                      metal_coords[1],
                                                      metal_coords[2]))
            # If we're within the distance expected, add this metal
            # to the list and move to the next
            if(distance - ionicRadiusDict[metal.type] <= radius):
                contacts.append(metal)
                break

    return contacts
