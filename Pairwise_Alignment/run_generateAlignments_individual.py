#!/usr/bin/env python3
import sys
from src.pairwiseGlobalAlignment_individual import pairwiseGlobalAlignment
import pandas as pd
import click
import datetime
import os
import multiprocessing
sys.path.append("src/")
default_outfile = "alignments_out_{date:%Y%m%d_%H%M%S}.csv"\
    .format(date=datetime.datetime.now())
default_outdir = "alignments_out_{date:%Y%m%d_%H%M%S}".format(
    date=datetime.datetime.now())


@click.command()
@click.option('--seqfile', default="21mers.csv",
              help="CSV file with 'sequences' and 'labels' columns")
@click.option('--matrixfile', default="matricies/HIJACK30",
              help="Text file with alignment matrix")
@click.option('--outfile', default=default_outfile,
              help="Name of the alignments outfile")
@click.option('--outdir', default=default_outdir,
              help="Name of the alignment outfolder")
@click.option('--ncores', default=multiprocessing.cpu_count())
def alignWrapper(seqfile, matrixfile, outfile, ncores, outdir):
    # Make output directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Import sequences
    sequences = pd.read_csv(seqfile)

    # Get the unique 21mer training data features and their associated labels
    kmerSequences = list(sequences["sequences"].unique())
    siteByStructureLabels = list(sequences["labels"].unique())

    # Align the sequences
    print("Running with {} cores".format(ncores))
    alignments = pairwiseGlobalAlignment(kmerSequences,
                                         siteByStructureLabels, outdir, ncores, matrixfile)

    # Set up a dataframe to view the pairwise global alignments
    alignmentsdf = pd.DataFrame(alignments[0])
    alignmentsdf.columns = siteByStructureLabels
    alignmentsdf["Site by Structure"] = siteByStructureLabels
    alignmentsdf.set_index("Site by Structure", inplace=True)

    alignmentsdf.to_csv(outfile)


if __name__ == '__main__':
    alignWrapper()
