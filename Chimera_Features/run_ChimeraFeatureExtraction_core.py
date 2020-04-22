#!/usr/bin/env chimera
# run_chimeraFeatureExtraction_core.py

"""
runChimeraFeatureExtraction_core.py

This script is a delagator. It takes in command line arguments, loop through each PDB
file passed by the user, and comandeers the main ChimeraFeatures package scripts
to get the features from that file, and then formats those features nicely for the user.

(This script and the associated package are formatted with black pypy.org/project/black)
"""

# Chimera-specific hackery below:
# Allows script to run without a chimera python interpreter, so users can run --help.
try:
    from src.chimeraFeatureExtraction import chimeraFeatureExtraction
    from chimera import runCommand as rc
except Exception as e:
    raise
    print(e)
    pass

import click
import datetime
import os
import sys
import pandas as pd
import time
from collections import deque
import traceback

# Set default names for output files, directories, etc.
default_outfile = "chimeraFeatureExtraction_out_{date:%Y%m%d_%H%M%S}.csv".format(
    date=datetime.datetime.now()
)
default_logfile = "chimeraFeatureExtraction_log_{date:%Y%m%d_%H%M%S}.txt".format(
    date=datetime.datetime.now()
)
default_outdir = "ChimeraOut_{date:%Y%m%d_%H%M%S}".format(date=datetime.datetime.now())
default_metalfile = "input/metalBindingSiteData_01-13-2017.txt"
default_sppider_binding = "input/sppider_binding.csv"
default_disopred_binding = "input/disopred_binding.csv"
default_disopred_disorder = "input/disopred_disorder.csv"


# Configure the command line parameters that this script can take in
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
    default="10",
    type=int,
    help="Number of re-attempts to generate surface computation after failure",
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
@click.option(
    "--outdir", default=default_outdir, help="Folder for all of the output files"
)
@click.option("--individual", default=True, help="Output individual files?")
@click.option("--threads", default=1, help="How many threads to run?")
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
    outdir,
    individual,
    threads,
):
    """ Uses command line arguments to loop through user-required pdb files.
    Then gets features for each pdb file passed and writes features to file.

    Arguments:
        pdbfile (str): Command line argument. Path to single pdb file.
        pdbdir (str): Command line argument. Path to directory of pdb files
        outfile (str): Command line argument. Path/name of output csv file
        radii (str): Space-delineated string of radii at which to compute
        bubble features.
    Returns:
        Void
    Effect:
        Writes a .csv file to "outfile" and (maybe) individual .csv files to outdir/individual
    """
    # Setup output directory (for .csv files)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if individual:
        if not os.path.exists(outdir + "/individual"):
            os.makedirs(outdir + "/individual")

    # Get the PDB files for this run, or throw an error if you can't find them
    # pdbdir overwrites pdbfile
    if pdbfile != "":
        pdb_files = [pdbfile]
    if pdbdir != "":
        pdb_files = [
            pdbdir + "/" + fn for fn in os.listdir(pdbdir) if fn.endswith(".pdb")
        ]
    if pdbfile == "" and pdbdir == "":
        raise ValueError("Must supply either a pdbfile or a pdbdir")

    # Just the main filename, no path
    pdb_files_short = [file.split("/")[-1] for file in pdb_files]

    print("Loading metal binding information into memory")
    # Get metal binding information, or default to empty if you cant find them
    try:
        metal_binding = get_metal_binding(metals_file)
    except Exception as e:
        metal_binding = {}
        log_message(
            logfile,
            "Metal binding file({}) not found: {}, default to empty dictionary".format(
                metals_file, e
            ),
        )

    # Get disopred disorder information, or default to empty if you can't find them
    # try:
    print("Loading disopred disorder information into memory")
    try:
        disopred_disorder = get_disopred_disorder(
            disopred_disorder_file, pdb_files_short
        )
    except Exception as e:
        disopred_binding = {}
        log_message(
            logfile,
            "DISOPRED binding file({}) not found: {}".format(disopred_binding_file, e),
        )

    # Get disopred binding information, or default to empty if you can't find them
    # try:
    print("Loading disopred binding information into memory")
    try:
        disopred_binding = get_disopred_binding(disopred_binding_file, pdb_files_short)
    except Exception as e:
        disopred_binding = {}
        log_message(
            logfile,
            "DISOPRED binding file({}) not found: {}".format(disopred_binding_file, e),
        )

    # Get sppider binding information, or default to empty if you can't find them
    # try:
    print("Loading sppider binding information into memory")
    try:
        sppider_binding = get_sppider_binding(sppider_binding_file, pdb_files_short)
    except Exception as e:
        sppider_binding = {}
        log_message(
            logfile,
            "SPPIDER binding file({}) not found: {}".format(sppider_binding_file, e),
        )

    log_message(logfile, "Beginning to process {} files\n".format(len(pdb_files)))

    # Note: Cannot figure out how to multithread the processing of these files. I think it's a
    #       bit weird because each requires a chimera instance, so multithreading would have to happen
    #       at a higher level, but then each higher level would use a lot of memory.

    # Now time to process all of the pdb files.
    pdb_files = deque(pdb_files)
    all_features = {}
    while len(pdb_files) >= 1:
        current_pdb = pdb_files.pop()
        try:
            log_message("chimera_progress.log", "Processing: {}\n".format(current_pdb))
            start_time = time.time()

            formatted_features = compute_features(
                current_pdb,
                metal_binding,
                disopred_disorder,
                disopred_binding,
                sppider_binding,
                attempts_limit,
                radii,
                logfile,
                individual,
                outdir,
            )
            all_features.update(formatted_features)
            with open("chimera_progress.log", "a") as progress_log:
                progress_log.write("Success processing: {}\n".format(current_pdb))
                progress_log.write(
                    "Time: %s seconds\n" % round(time.time() - start_time, 2)
                )

            log_message(logfile, "Feature extraction success: {}\n".format(current_pdb))

        except Exception as e:
            # Write error to logs
            log_message(
                logfile,
                "Feature extraction failed {},{},{},{}\n\n".format(
                    current_pdb, e, sys.exc_info(), traceback.format_exc()
                ),
            )
            log_message("chimera_progress.log", "Failure processing: {}\n".format(current_pdb))
            log_message("chimera_errors.log", "{}\t{}\n".format(current_pdb, traceback.format_exc().strip()))

            continue

    # We really just need features from RPKT atoms, so put those into one dictionary
    # and then write them to file
    rpkt_features = {}
    for atom_name, features in all_features.items():
        if atom_name.split()[1] in ["R", "P", "K", "T"]:
            rpkt_features[atom_name] = features

    write_all_features(rpkt_features, outdir + "/" + outfile)


def compute_features(
    file,
    metal_binding,
    disopred_disorder,
    disopred_binding,
    sppider_binding,
    attempts_limit,
    radii,
    logfile,
    individual,
    outdir,
):
    try:
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

        # Extract features using the chimeraFeatureExtraction script
        features = chimeraFeatureExtraction(
            file,
            [int(radius) for radius in radii.split()],
            metal_binding_sites,
            disopred_disorder_map,
            disopred_binding_map,
            sppider_binding_map,
            logfile,
            attempts_limit,
        )

        # Re-format features so that they're easier to write
        formatted_features = format_single_features(features, file)

        if individual:
            # Write the formatted features to their indiv files if parameter set
            single_rpkt_features = {}
            for atom_name, features in formatted_features.items():
                if atom_name.split()[1] in ["R", "P", "K", "T"]:
                    single_rpkt_features[atom_name] = features

            write_all_features(
                single_rpkt_features,
                "{}/individual/{}.csv".format(
                    outdir, file_name_short.replace(".pdb", "")
                ),
            )

        rc("close all")

        return formatted_features

    # Catch exception so we can close all windows before throwing the error
    except Exception as e:
        rc("close all")
        raise


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


def log_message(logfile, message):
    with open(logfile, "a") as log:
        log.write(message)


if __name__ == "__main__":
    featureWrapper()
