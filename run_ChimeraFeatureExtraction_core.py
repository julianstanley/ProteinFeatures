"""
runChimeraFeatureExtraction_core.py

This script is meant to take in command line arguments and use those to extract features
from .pdb protein files.

This script and the associated package are formatted with black pypy.org/project/black
"""

# Use a try block for importing chimeraFeatureExtraction so that you
# can run ./run_ChimeraFeatureExtraction_core.py --help
# (Since Chimera can't be imported without running via chimera)
try:
    from src.chimeraFeatureExtraction import chimeraFeatureExtraction
    from chimera import runCommand as rc
except Exception as e:
    print(e)
    pass

# Use the click module to set command line arguments
import click
import datetime
import os
import pandas as pd
import time

# Set default names for the log and outfile (exported) as well as the
# metal binding site file (imported)
default_outfile = "chimeraFeatureExtraction_out_{date:%Y%m%d_%H%M%S}.csv".format(
    date=datetime.datetime.now()
)
default_logfile = "chimeraFeatureExtraction_log_{date:%Y%m%d_%H%M%S}.txt".format(
    date=datetime.datetime.now()
)
default_metalfile = "input/metalBindingSiteData_01-13-2017.txt"
default_sppider_binding = "input/sppider_binding.csv"
default_disopred_binding = "input/disopred_binding.csv"
default_disopred_disorder = "input/disopred_disorder.csv"


# Set command-line parameters, to be passed to the "featureWrapper" function
@click.command()
@click.option("--pdbfile", default="", help="(Optional) PDB file to analyze")
@click.option("--pdbdir", default="", help="(Optional) directory of PDB files")
@click.option(
    "--outfile",
    default=default_outfile,
    help="(Optional) filename to save .csv output file",
)
@click.option(
    "--radii",
    default="5 8 10 12",
    help="""Space seperated ints representing radii
               at which to take bubble features""",
)
@click.option(
    "--attempts_limit",
    default="1",
    help="Number of re-attempts to generate features after failure",
)
@click.option(
    "--logfile", default=default_logfile, help="The log file for errors/comments/etc."
)
@click.option(
    "--metals_file",
    default=default_metalfile,
    help="Path containing the metal binding data file",
)
@click.option(
    "--disopred_disorder_file",
    default=default_disopred_disorder,
    help="Path containing the disopred disorder data file",
)
@click.option(
    "--disopred_binding_file",
    default=default_disopred_binding,
    help="Path containing the disopred protein binding data file",
)
@click.option(
    "--sppider_binding_file",
    default=default_sppider_binding,
    help="Path containing the SPPIDER protein binding data file",
)
def featureWrapper(
    pdbfile,
    pdbdir,
    outfile,
    radii,
    attempts_limit,
    logfile,
    metals_file,
    disopred_disorder_file,
    disopred_binding_file,
    sppider_binding_file,
):
    """ Takes in the command line arguments, checks to make sure pdb files
    were passed, and then loops through each pdb file, grabs features from
    chimeraFeatureExtraction.py, formats those features with
    format_single_features, and then writes those features
    with write_all_features()
    Arguments:
        pdbfile (str): Command line argument. Path to single pdb file.
        pdbdir (str): Command line argument. Path to directory of pdb files
        outfile (str): Command line argument. Path/name of output csv file
        radii (str): Space-delineated string of radii at which to compute
        bubble features.
    Returns:
        Void
    Effect:
        Writes a .csv file to "outfile"
    """
    # If we fail, how many times are we going to retry? This is passed as a string
    # from the command line, but should be an int.
    attempts_limit = int(attempts_limit)

    # Get the PDB files for this run, or throw an error if you can't find them
    if pdbfile != "":
        pdb_files = [pdbfile]
    if pdbdir != "":
        pdb_files = [
            pdbdir + "/" + fn for fn in os.listdir(pdbdir) if fn.endswith(".pdb")
        ]
    if pdbfile == "" and pdbdir == "":
        raise ValueError("Must supply either a pdbfile or a pdbdir")

    pdb_files_short = [file.split("/")[-1] for file in pdb_files]

    print("Loading metal binding information into memory")

    # Get metal binding information, or default to empty if you cant find them
    try:
        metal_binding = get_metal_binding(metals_file)
    except Exception as e:
        metal_binding = {}
        with open(logfile, "w+") as log:
            log.write("Metal binding file({}) not found: {}".format(metals_file, e))

    all_features = {}

    # Get disopred disorder information, or default to empty if you can't find them
    # try:
    print("Loading disopred disorder information into memory")
    disopred_disorder = get_disopred_disorder(disopred_disorder_file, pdb_files_short)
    # except Exception as e:
    #     disopred_disorder = {}
    #     with open(logfile, "w+") as log:
    #         log.write(
    #             "DISOPRED disorder file({}) not found: {}".format(
    #                 disopred_disorder_file, e
    #             )
    #         )

    # Get disopred binding information, or default to empty if you can't find them
    # try:
    print("Loading disopred binding information into memory")
    disopred_binding = get_disopred_binding(disopred_binding_file, pdb_files_short)
    # except Exception as e:
    #     disopred_binding = {}
    #     with open(logfile, "w+") as log:
    #         log.write(
    #             "DISOPRED binding file({}) not found: {}".format(
    #                 disopred_binding_file, e
    #             )
    #         )

    # Get sppider binding information, or default to empty if you can't find them
    # try:
    print("Loading sppider binding information into memory")
    sppider_binding = get_sppider_binding(sppider_binding_file, pdb_files_short)
    # except Exception as e:
    #     sppider_binding = {}
    #     with open(logfile, "w+") as log:
    #         log.write(
    #             "SPPIDER binding file({}) not found: {}"
    # .format(sppider_binding_file, e)
    #         )

    with open("chimera_progress.log", "a") as progress_log:
        progress_log.write("Beginning to process {} files\n".format(len(pdb_files)))

    for file in pdb_files:
        with open("chimera_progress.log", "a") as progress_log:
            progress_log.write("Processing: {}\n".format(file))
        start_time = time.time()

        file_name_short = file.split("/")[-1]
        # Get the metal binding sites associated with the file, if they exist
        if file_name_short in metal_binding:
            metal_binding_sites = metal_binding[file_name_short]
        else:
            metal_binding_sites = []

        # Get the disopred disorder info associated with the file, if it exists
        if file_name_short in disopred_disorder:
            disopred_disorder_map = disopred_disorder[file_name_short]
        else:
            disopred_disorder_map = []

        # Get the disopred binding info associated with the file, if it exists
        if file_name_short in disopred_binding:
            disopred_binding_map = disopred_binding[file_name_short]
        else:
            disopred_binding_map = []

        # Get the sppider binding info associated with the file, if it exists
        if file_name_short in sppider_binding:
            sppider_binding_map = sppider_binding[file_name_short]
        else:
            sppider_binding_map = []

        attempts = 0
        while attempts < attempts_limit:
            try:
                # Extract features using the chimeraFeatureExtraction script
                features = chimeraFeatureExtraction(
                    file,
                    [int(radius) for radius in radii.split()],
                    metal_binding_sites,
                    disopred_disorder_map,
                    disopred_binding_map,
                    sppider_binding_map,
                    logfile,
                )

                # Re-format features so that they're easier to write
                all_features.update(format_single_features(features, file))

                # We got the features we were looking for, so no need to repeat
                with open("chimera_progress.log", "a") as progress_log:
                    progress_log.write("Success processing: {}\n".format(file))
                    progress_log.write(
                        "Time: %s seconds\n" % round(time.time() - start_time, 2)
                    )
                rc("close all")
                break
            except Exception as e:
                with open(logfile, "w+") as log:
                    log.write("Feature extraction failed {}, {}\n".format(file, e))
                with open("chimera_progress.log", "a") as progress_log:
                    progress_log.write("Failure processing: {}\n".format(file))
                rc("close all")
                break

    # We really just need features from RPKT atoms, so put those into one dictionary
    # and then write them to file
    rpkt_features = {}
    for atom_name, features in all_features.items():
        if atom_name.split()[1] in ["R", "P", "K", "T"]:
            rpkt_features[atom_name] = features

    write_all_features(rpkt_features, outfile)


def get_metal_binding(metal_file_path):
    """ Parses a file with metal binding data
    Arguments:
        metal_file_path (str): A directory path to the file with metal binding data
    Returns:
        file_to_metals (dict): A dictionary of metal binding data
            Format: {File Name: str, Site Number: str, Metal ID: str, Residues: str}
    Effect: None
    """
    file_to_metals = {}
    with open(metal_file_path, "r") as f:
        for line in f:
            # The metal binding file contains:
            # Protein file name, site number, metal (Zn, Mg, etc)
            # and then residue (D127, Y125, etc.)
            file_name, site_num, metals, residues = line.split()
            if file_name not in file_to_metals:
                file_to_metals[file_name] = []

            # Entries should have a unique metal identifier
            for metal_id in metals.split(","):
                file_to_metals[file_name].append(
                    {
                        "File Name": file_name,
                        "Site Number": site_num,
                        "Metal ID": metal_id,
                        "Residues": residues,
                    }
                )

    return file_to_metals


def get_disopred_disorder(disopred_disorder_file, pdb_files_short):
    file_to_disorder = {}
    with open(disopred_disorder_file, "r") as f:
        for line in f:
            line = line.strip()
            # The disopred disorder file contains:
            # Protein file name, amino acid, position, disorder call, disorder score
            file_name, aa, pos, disorder_call, disorder_score = line.split(",")
            if file_name in pdb_files_short:
                if file_name not in file_to_disorder:
                    file_to_disorder[file_name] = []

                try:
                    file_to_disorder[file_name].append(
                        {
                            "File Name": file_name,
                            "aa": aa,
                            "pos": pos,
                            "Disorder Call": int(disorder_call),
                            "Disorder Score": float(disorder_score),
                        }
                    )
                except Exception as e:
                    continue

    return file_to_disorder


def get_disopred_binding(disopred_binding_file, pdb_files_short):
    file_to_binding = {}
    with open(disopred_binding_file, "r") as f:
        for line in f:
            # The disopred disorder file contains:
            # Protein file name, amino acid, position, binding call, binding score
            file_name, aa, pos, binding_call, binding_score = line.split(",")
            if file_name in pdb_files_short:
                if file_name not in file_to_binding:
                    file_to_binding[file_name] = []

                try:
                    file_to_binding[file_name].append(
                        {
                            "File Name": file_name,
                            "aa": aa,
                            "pos": pos,
                            "Binding Call": int(binding_call),
                            "Binding Score": float(binding_score),
                        }
                    )
                except Exception as e:
                    continue

    return file_to_binding


def get_sppider_binding(sppider_binding_file, pdb_files_short):
    file_to_binding = {}
    with open(sppider_binding_file, "r") as f:
        for line in f:
            # The disopred disorder file contains:
            # Protein file name, amino acid, position, binding call
            file_name, aa, pos, binding_call = line.split(",")
            if file_name in pdb_files_short:
                if file_name not in file_to_binding:
                    file_to_binding[file_name] = []

                try:
                    file_to_binding[file_name].append(
                        {
                            "File Name": file_name,
                            "aa": aa,
                            "pos": pos,
                            "Binding Call": int(binding_call),
                        }
                    )
                except Exception as e:
                    continue

    return file_to_binding


def format_single_features(single_features, file_name):
    """ Puts global, bubble, and residue-specific features into a single,
    uniformally-formatted dictionary.
    Arguments:
        single_features(tuple, length 3): The global, residue-level, and bubble-level
        features of a protein
        file_name(str): The name of the protein file name that we're processing
    Returns:
        atom_to_features(dict): A dictionary mapping atoms to their features
    """
    global_features = single_features[0]
    atom_features = single_features[1]
    atom_radius_features = single_features[2]

    atom_to_features = {}

    sum_excluded = [
        "circularVariance",
        "sumMetals",
        "CA",
        "CO",
        "CU",
        "FE",
        "K",
        "MG",
        "MN",
        "MO",
        "NA",
        "NI",
        "ZN",
    ]

    for atom, radius_to_bubble_features in atom_radius_features.items():
        features = {}
        features["Protein"] = file_name.split("/")[-1]
        for feature_name, feature_value in global_features.items():
            if (
                all([excl not in feature_name for excl in sum_excluded])
                or feature_name == "RKPT"
            ):
                if feature_name == "contacts":
                    features["residues_global_sum"] = feature_value
                else:
                    features["{}_global_sum".format(feature_name)] = feature_value
            else:
                features["{}_global".format(feature_name)] = feature_value

        for feature_name, feature_value in atom_features[atom].items():
            features[feature_name] = feature_value

        for radius, bubble_features in radius_to_bubble_features.items():
            for bubble_feature_name, feature_value in bubble_features.items():
                if (
                    all([excl not in bubble_feature_name for excl in sum_excluded])
                    or bubble_feature_name == "RKPT"
                ):
                    feature_name = "{bubble}_{radius}A_sum".format(
                        bubble=bubble_feature_name, radius=radius
                    )
                else:
                    feature_name = "{bubble}_{radius}A".format(
                        bubble=bubble_feature_name, radius=radius
                    )
                features[feature_name] = feature_value
                # if 'contacts' not in feature_name:
                #     features[feature_name.replace("sum", "avg")] =\
                #         feature_value / bubble_features["contacts"]

        atom_to_features[atom] = features

    return atom_to_features


def write_all_features(all_features, outfile):
    """ Writes Chimera features to the given outfile .csv
    Arguments:
        all_features (dictionary)
            format: {Atom Name --> Feature Object (NamedTuple),
            'Protein' --> Protein File Name (Str)}: All Chimera features
        outfile (string) A filename for the features .csv file
    Returns:
        Void
    Effect:
        Creates a file at 'outfile'
    """
    features_df = pd.DataFrame(all_features).transpose()
    features_df.to_csv(outfile)


if __name__ == "__main__":
    featureWrapper()
