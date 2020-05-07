#!/usr/bin/env python3
# Description: a relatively hastily-written python script to format the output files from the chimera pipeline.
# Please feel free to further modularize and clean things up! :) Or, if you want me to clean it up,
# email me at julianstanleyA@gmail.com
# Call `--help` for more info.
import click
import pandas as pd
import datetime

# Set default names for output files, directories, etc.
default_outfile = "chimeraFormat_out_{date:%Y%m%d_%H%M%S}.csv".format(
    date=datetime.datetime.now()
)
default_columnsfile = "desired_columns.txt"

# Configure the command line parameters that this script can take in
@click.command()
@click.option("--features",
              default="",
              help="Path to the unformatted features file. Required.")
@click.option("--columns",
              default=default_columnsfile,
              help="A file with the columns to extact/format, one per line (Default: desired_columns.txt)")
@click.option("--outfile",
              default=default_outfile,
              help="The output file (default: chimeraFormat_out_[date].csv)")
def main(features, columns, outfile):
    # Load the unformatted data
    unformatted_features = pd.read_csv(features, low_memory=False)
    # Remove any unnamed columns from unformatted data
    unformatted_features = unformatted_features[[
        x for x in unformatted_features.columns if "Unnamed" not in str(x)]]

    # Load the desired columns
    with open(columns, "r") as f:
        desired_columns = [x.strip() for x in f]

    # Add a "Site by Structure" column to the unformatted data
    position = []
    for x in unformatted_features["Position"]:
        try:
            position.append(x.split(".")[0])
        except Exception as e:
            position.append(str(x))

    aa = unformatted_features["AA"]
    protein = unformatted_features["Protein"]

    unformatted_features["Site by Structure"] = protein + "_" + aa + position

    # Does this dataset include a "None" column?
    nones = any(["_None_" in str(x) for x in unformatted_features.columns])

    # Add "averages" feature columns to data
    for feature in [y for y in desired_columns if y not in unformatted_features and y.split("_")[-1] == "avg"]:
        associated_sum = feature.replace("avg", "sum")
        sum_new = unformatted_features[associated_sum]
        if "global" in feature:
            if nones:
                associated_none = associated_sum.replace(
                    "_global_", "_None_global_")
                if associated_none in unformatted_features.columns:
                    res = unformatted_features["residues_global_sum"] - \
                        unformatted_features[associated_none]
                else:
                    associated_none = associated_sum.replace(
                        "_global_sum", "_None_global")
                    res = unformatted_features["residues_global_sum"] - \
                        unformatted_features[associated_none]
            else:
                res = unformatted_features["residues_global_sum"]
        else:
            distance = feature.split("_")[1]
            if nones:
                res = unformatted_features[f'contacts_{distance}_sum'] + \
                    1 - unformatted_features[f'contacts_None_{distance}_sum']
            else:
                res = unformatted_features[f'contacts_{distance}_sum'] + 1

        average_new = sum_new / res
        unformatted_features[feature] = average_new

    formatted_features = pd.DataFrame()
    formatted_features["Site by Structure"] = unformatted_features["Site by Structure"]
    for colname in list(desired_columns):
        try:
            formatted_features[colname] = unformatted_features[colname]
        except Exception as e:
            print(f"Can't find column: {colname}, so leaving out. Error: {e}")
            continue

    formatted_features.to_csv(outfile, index=None)


if __name__ == "__main__":
    main()
