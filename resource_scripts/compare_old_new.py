#!/usr/bin/env python3
import pandas as pd

# Grab the new and old feature frames
old = pd.read_csv("../data/old_ref.csv")
new = pd.read_csv("../new_compare.csv")

# logging
log = open("logfile.txt", "w")

# Add a Site by Structure column to the new dataframe
new["Site by Structure"] = new.apply(
    lambda x: "{protein}_{res}{loc}".format(
        protein=x["Protein"], res=x["AA"], loc=x["Position"].split(".")[0]
    ),
    axis=1,
)

outfile = open("old_new_compare.csv", "w")
outfile.write("Feature,Old,New,Match\n")

for index, row in old.iterrows():
    try:
        new_match = new[new["Site by Structure"] == row["Site by Structure"]]
        if not new_match.empty:
            outfile.write(
                "New Structure,{},{},TRUE\n".format(
                    row["Site by Structure"], row["Site by Structure"]
                )
            )
            for old_name in list(row.keys()):
                if old_name in list(new_match):
                    feature = old_name
                    old_feature = row[old_name]
                    new_feature = new_match[old_name].values[0]
                    if isinstance(new_feature, (int, float)) and not isinstance(
                        new_feature, bool
                    ):
                        match = (
                            abs(float(old_feature) - float(new_feature)) < 0.001
                        )
                    else:
                        match = old_feature == new_feature
                    outfile.write(f"{feature},{old_feature},{new_feature},{match}\n")

                else:
                    log.write("Not matched: {}\n".format(old_name))
    except Exception as e:
        print(e)
outfile.close()
