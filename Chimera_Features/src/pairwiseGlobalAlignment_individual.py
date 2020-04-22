from Bio import pairwise2
import pandas as pd
from scipy.sparse import csr_matrix
import csv
import multiprocessing
from functools import partial


def pairwiseGlobalAlignment(seqs, siteByStructureLabels, outdir, ncores, matDir="../matricies/HIJACK30"):
    """ Align a set of amino acid sequences pairwise/amongst themselves

    Args:
        seqs (List[str]): A list of amino acid sequences of equal length
        for pairwise alignment

        matDir (str, optional): The directory containing the alignment matrix

    Returns:
        List[ndarray(float)]: Two returns in a list
        [0]: Alignment scores (ndarray(float)) of size len(seqs) x len(seqs).
        The i,j entry of this array represents the alignment between seqs[i]
        and seqs[j].
        [1]: Lengths of size len(seqs) x len(seqs).
        The i,j entry of this array represents the minimum length of seqs[i]
        and seqs[j].

    Notes:
    This uses the HIJACK30 matrix, which must be included in the working directory
    of this function. The matrix is in standard NCBI format and is a modified BLOSUM30
    matrix such that a '*' to '*' match has a score of 1000 and all matches with 'X'
    have a score of 0.

    """

    # Read the HIJACK30 matrix from file, then collapse into a dictionary format
    HIJACK30_df = pd.read_csv(matDir, sep="\t")
    HIJACK30 = {}
    for columnName in HIJACK30_df.columns:
        for rowIndex in range(24):
            HIJACK30[(columnName, HIJACK30_df.columns[rowIndex])] = HIJACK30_df[
                columnName
            ][rowIndex]

    # Store the number of sequences
    numSeqs = len(seqs)

    # Initalize return results
    scores = csr_matrix((numSeqs, numSeqs))
    lengths = csr_matrix((numSeqs, numSeqs))

    # Set column labels
    column_labels = siteByStructureLabels[:]
    column_labels.insert(0, "Site By Structure")

    # Compare all sequences with themselves (without replacement,
    # so like an upper triangular matrix)
    pool = multiprocessing.Pool(ncores)
    pool.map(partial(get_alignment_row, numSeqs=numSeqs,
                     seqs=seqs, HIJACK30=HIJACK30,
                     siteByStructureLabels=siteByStructureLabels,
                     outdir=outdir, column_labels=column_labels), range(numSeqs))

    return [scores, lengths]


def get_alignment_row(i, numSeqs, seqs, HIJACK30, siteByStructureLabels, outdir,
                      column_labels):
    row_label = siteByStructureLabels[i]
    alignment_row = [float('NaN')] * numSeqs
    lengths = [float('NaN')] * numSeqs
    for j in range(numSeqs):
        if i <= j:
            seq1 = list(seqs[i])
            seq1Len = len(seq1)
            seq2 = list(seqs[j])
            seq2Len = len(seq2)

            # Remember the minimum length of the two sequences
            lengths[j] = min(seq1Len, seq2Len)

            # Remember central residues in the two sequences and their indicies
            seq1CenterIndex = seq1Len // 2
            seq2CenterIndex = seq2Len // 2
            seq1Center = seq1[seq1CenterIndex]
            seq2Center = seq2[seq2CenterIndex]

            # Replace central residues with '*' to force residue alignment
            seq1[seq1CenterIndex] = "*"
            seq2[seq2CenterIndex] = "*"

            # Calculate the alignment score for both the sequences and the center
            # residue using the HIJACK30 matrix.
            # The gap opening and extension penalities of -8
            # and penalize_end_gaps=False parameter
            # let the pairwise.align function with similarly (within +/-1.5E-14)
            # to the glocal parameter in MATLAB on 979 unique 21-mers
            score = pairwise2.align.globalds(
                "".join(seq1),
                "".join(seq2),
                HIJACK30,
                -8,
                -8,
                penalize_end_gaps=False,
                score_only=True,
            )
            scoreCenter = pairwise2.align.globalds(
                seq1Center,
                seq2Center,
                HIJACK30,
                -8,
                -8,
                penalize_end_gaps=False,
                score_only=True,
            )

            # The alignment scores should be scaled by 0.2
            scoreScaled = 0.2 * score
            scoreScaledCenter = 0.2 * scoreCenter

            # Correct the alignment score (remove the *-* alignment and add a
            # center-center alignment)
            scoreCorrected = scoreScaled + scoreScaledCenter - 1000 * 0.2
            alignment_row[j] = scoreCorrected

    outfile = "{}/{}.csv".format(outdir, row_label)

    with open(outfile, 'w') as f:
        alignment_row.insert(0, row_label)
        writer = csv.writer(f)
        writer.writerow(column_labels)
        writer.writerow(alignment_row)

    return alignment_row
