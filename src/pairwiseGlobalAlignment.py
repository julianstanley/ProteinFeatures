import numpy as np
from Bio import pairwise2
import pandas as pd
from itertools import combinations_with_replacement


def pairwiseGlobalAlignment(seqs, matDir="../matricies/HIJACK30"):
    """ Align a set of amino acid sequences pairwise/amongst themselves

    Args:
        seqs (List[str]): A list of amino acid sequences of equal length 
        for pairwise alignment
        
        matDir (str, optional): The directory containing the alignment matrix

    Returns:
        List[ndarray(float)]: Two returns in a list
        [0]: Alignment scores (ndarray(float)) of size len(seqs) x len(seqs). 
        The i,j entry of this array represents the alignment between seqs[i] and seqs[j].
        [1]: Lengths of size len(seqs) x len(seqs). 
        The i,j entry of this array represents the minimum length of seqs[i] and seqs[j].
    
    Notes:
    This uses the HIJACK30 matrix, which must be included in the working directory of this function. 
    The matrix is in standard NCBI format and is a modified BLOSUM30 matrix such that a '*' to '*' 
    match has a score of 1000 and all matches with 'X' have a score of 0. 
    
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
    scores = np.zeros((numSeqs, numSeqs))
    lengths = np.zeros((numSeqs, numSeqs))

    # Compare all sequences with themselves (without replacement, so like an upper triangular matrix)
    for i, j in combinations_with_replacement(range(numSeqs), 2):
        # Get the two sequences to align
        # Note: working with strings as lists, since lists are mutable
        seq1 = list(seqs[i])
        seq1Len = len(seq1)
        seq2 = list(seqs[j])
        seq2Len = len(seq2)

        # Remember the minimum length of the two sequences
        lengths[i, j] = min(seq1Len, seq2Len)
        lengths[j, i] = min(seq1Len, seq2Len)

        # Remember central residues in the two sequences and their indicies
        seq1CenterIndex = seq1Len // 2
        seq2CenterIndex = seq2Len // 2
        seq1Center = seq1[seq1CenterIndex]
        seq2Center = seq2[seq2CenterIndex]

        # Replace central residues with '*' to force residue alignment
        seq1[seq1CenterIndex] = "*"
        seq2[seq2CenterIndex] = "*"

        # Calculate the alignment score for both the sequences and the center residue using the HIJACK30 matrix
        # The gap opening and extension penalities of -8 and penalize_end_gaps=False parameter let the pairwise.align
        # function with similarly (within +/-1.5E-14) to the glocal parameter in MATLAB on 979 unique 21-mers
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

        # Correct the alignment score (remove the *-* alignment and add a center-center alignment)
        scoreCorrected = scoreScaled + scoreScaledCenter - 1000 * 0.2
        scores[i, j] = scoreCorrected
        scores[j, i] = scoreCorrected

    return [scores, lengths]
