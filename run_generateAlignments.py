#!/usr/bin/env python3
import sys
from src.pairwiseGlobalAlignment import pairwiseGlobalAlignment
import pandas as pd
import click
import datetime
sys.path.append("src/")
default_outfile = "alignments_out_{date:%Y%m%d_%H%M%S}.csv"\
    .format(date=datetime.datetime.now())


@click.command()
@click.option('--seqfile', default="21mers.csv",
              help="CSV file with 'sequences' and 'labels' columns")
@click.option('--matrixfile', default="matricies/HIJACK30",
              help="Text file with alignment matrix")
@click.option('--outfile', default=default_outfile,
              help="Name of the alignments outfile")
def alignWrapper(seqfile, matrixfile, outfile):
    # Import sequences
    sequences = pd.read_csv(seqfile)

    # Get the unique 21mer training data features and their associated labels
    kmerSequences = list(sequences["sequences"].unique())
    siteByStructureLabels = list(sequences["labels"].unique())

    # Align the sequences
    alignments = pairwiseGlobalAlignment(kmerSequences, matrixfile)

    # Set up a dataframe to view the pairwise global alignments
    alignmentsdf = pd.DataFrame(alignments[0])
    alignmentsdf.columns = siteByStructureLabels
    alignmentsdf["Site by Structure"] = siteByStructureLabels
    alignmentsdf.set_index("Site by Structure", inplace=True)

    alignmentsdf.to_csv(outfile)


if __name__ == '__main__':
    alignWrapper()
