# Use a try block for importing chimeraFeatureExtraction so that you
# can run ./run_ChimeraFeatureExtraction_core.py --help
# (Since Chimera can't be imported without running via chimera)
try:
    from src.chimeraFeatureExtraction import chimeraFeatureExtraction
    from chimera import runCommand as rc
    from src.ExtendedChimera.ChimeraUtils import *
except Exception as e:
    print(e)
    pass
# Use the click module to set command line arguments
import click
import datetime
import os
import pandas as pd
default_outfile = "chimeraFeatureExtraction_out_{date:%Y%m%d_%H%M%S}.csv"\
    .format(date=datetime.datetime.now())
default_logfile = "chimeraFeatureExtraction_log_{date:%Y%m%d_%H%M%S}.txt"\
    .format(date=datetime.datetime.now())
default_metalfile = "input/metalBindingSiteData_01-13-2017.txt"


@click.command()
@click.option("--pdbfile", default="", help="(Optional) PDB file to analyze")
@click.option("--pdbdir", default="", help="(Optional) directory of PDB files")
@click.option("--outfile", default=default_outfile,
              help="(Optional) filename to save .csv output file")
@click.option("--radii", default="5 8 10 12",
              help="""Space seperated ints representing radii
               at which to take bubble features""")
@click.option("--attempts_limit", default="5",
              help="Number of re-attempts to generate features after failure")
@click.option("--logfile", default=default_logfile,
              help="The log file for errors/comments/etc.")
@click.option("--metals_file", default=default_metalfile,
              help="Path containing the metal binding data file")
def featureWrapper(pdbfile, pdbdir, outfile, radii, attempts_limit, logfile,
                   metals_file):
    ''' Takes in the command line arguments, checks to make sure pdb files
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
    '''
    attempts_limit = int(attempts_limit)

    if pdbfile != "":
        pdb_files = [pdbfile]
    if pdbdir != "":
        pdb_files = [pdbdir + "/" + fn
                     for fn in os.listdir(pdbdir) if fn.endswith(".pdb")]
    if (pdbfile == "" and pdbdir == ""):
        raise ValueError("Must supply either a pdbfile or a pdbdir")

    # Get metal binding information, or empty dict if not found
    try:
        metal_binding = get_metal_binding(metals_file)
    except Exception as e:
        metal_binding = {}
        with open(logfile, 'w+') as log:
            log.write(
                "Metal binding file ({}) not found: {}".format(metals_file, e))

    all_features = {}
    for file in pdb_files:
        file_name_short = file.split("/")[-1]
        # Get the metal binding sites associated with the file, if it exists
        if file_name_short in metal_binding:
            metal_binding_sites = metal_binding[file_name_short]
        else:
            metal_binding_sites = []

        # The attempts limit was passed via command line as a string
        # attempts_limit = int(attempts_limit)
        attempts = 0
        while attempts < attempts_limit:
            try:
                features = chimeraFeatureExtraction(
                    file, [int(radius) for radius in radii.split()],
                    metal_binding_sites)

                all_features.update(format_single_features(features, file))
                break

            except Exception as e:
                rc("close all")
                with open(logfile, 'w+') as log:
                    log.write("Begin exception log:")
                    save_and_clear_reply_log(logfile)
                    log.write("""Feature extraction for {} failed {} times.
                    The most recent exception thrown was {}""".format(
                        file, attempts_limit, e))
                # continue
                break

    rpkt_features = {}
    for atom_name, features in all_features.items():
        if atom_name.split()[1] in ["R", "P", "K", "T"]:
            rpkt_features[atom_name] = features

    write_all_features(rpkt_features, outfile)


def get_metal_binding(metal_site_file):
    '''
    Returns: TODO: {File --> [List-of-Dicts]}
    Where the list of dicts correspond to each metal binding site
    in that file
    Where [of-Dicts] is in format: {Metal --> 'Zn, Ca, etc',
    Residues --> 'C208, H1, etc'}
    '''
    file_to_metals = {}
    with open(metal_site_file, 'r') as f:
        for line in f:
            # The metal binding file contains:
            # Protein file name, site number, metal (Zn, Mg, etc)
            # and then residue (D127, Y125, etc.)
            file_name, site_num, metals, residues = line.split()
            if file_name not in file_to_metals:
                file_to_metals[file_name] = []

            # Entries should have a unique metal identifier
            for metal_id in metals.split(","):
                file_to_metals[file_name].append({'File Name': file_name,
                                                  'Site Number': site_num,
                                                  'Metal ID': metal_id,
                                                  'Residues': residues
                                                  })

    return file_to_metals


def write_all_features(all_features, outfile):
    ''' Writes Chimera features to the given outfile .csv
    Arguments:
        all_features (dictionary)
            format: {Atom Name --> Feature Object (NamedTuple),
            'Protein' --> Protein File Name (Str)}: All Chimera features
        outfile (string) A filename for the features .csv file
    Returns:
        Void
    Effect:
        Creates a file at 'outfile'
    '''
    features_df = pd.DataFrame(all_features).transpose()
    features_df.to_csv(outfile)


def format_single_features(single_features, file_name):
    ''' Puts global, bubble, and residue-specific features into a single,
    uniformally-formatted dictionary.
    Arguments:
    Returns:
    '''
    global_features = single_features[0]
    residue_features = single_features[1]
    atom_radius_features = single_features[2]

    atom_to_features = {}

    for atom, radius_to_bubble_features\
            in atom_radius_features.items():
        features = {}
        features['Protein'] = file_name.split("/")[-1]
        for feature_name, feature_value in global_features.items():
            features["{}_global_sum".format(feature_name)] = feature_value
            if 'contacts' not in feature_name:
                features["{}_global_avg".format(feature_name)] = \
                    feature_value / global_features["contacts"]

        for feature_name, feature_value in residue_features[atom].items():
            features[feature_name] = feature_value

        for radius, bubble_features in radius_to_bubble_features.items():
            for bubble_feature_name, feature_value \
                    in bubble_features.items():
                feature_name = "{bubble}_{radius}A_sum".format(
                    bubble=bubble_feature_name, radius=radius)
                features[feature_name] = feature_value
                if 'contacts' not in feature_name:
                    features[feature_name.replace("sum", "avg")] =\
                        feature_value / bubble_features["contacts"]
        atom_to_features[atom] = features

    return(atom_to_features)


if __name__ == "__main__":
    featureWrapper()
