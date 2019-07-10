import numpy as np
from Bio.SubsMat import MatrixInfo as matlist
import nwalign3 as nw
import math

def pairwiseGlobalAlignment(seqs, matDir = ""):
    ''' Align a set of amino acid sequences pairwise/amongst themselves

    Args:
        seqs (List[str]): A list of amino acid sequences of equal length 
        for pairwise alignment
        
        matDir (str, optional): The directory containing the alignment matrix

    Returns:
        ndarray(float): A matrix of alignment 
    
    Notes:
    This uses the HIJACK30 matrix, which must be included 
    in the working directory of this function. The matrix
    is in standard NCBI format and is a modified BLOSUM30
    matrix such that a '*' to '*' match has a score of 1000
    and all matches with 'X' have a score of 0. 
    
    '''
    
    # Store the number of sequences
    numSeqs = len(seqs)
    
    # Initalize return results
    scores = np.zeros((numSeqs, numSeqs))
    lengths = np.zeros((numSeqs, numSeqs))
            
    # Note: using a non-parallelized for loop. If necessary, use 
    # multiprocessing module to speed up. 
    # Nested loop to compare all seqs against eachother
    # Exclude i=j to not compare sequence with itself
    for i in range(numSeqs):
        for j in range(numSeqs):
            # Get the two sequences to align
            # (and their lengths for references
            # Note: work with strings as mutable lists
            seq1 = list(seqs[i])
            seq1Len = len(seq1)
            seq2 = list(seqs[j])
            seq2Len = len(seq2)
            
            # Report back the minimum length of the two sequences
            lengths[i,j] = min(seq1Len, seq2Len)

            # Remember the number of the central residues in the two sequences
            seq1CenterIndex = math.ceil(seq1Len/2)
            seq2CenterIndex = math.ceil(seq2Len/2)
            
            # Parse central residues:
            seq1Center = seq1[seq1CenterIndex]
            seq2Center = seq2[seq2CenterIndex]
            
            # Replace central residues with '*' to force residue alignment
            seq1[math.ceil(seq1Len/2)] = '*'
            seq2[math.ceil(seq2Len/2)] = '*'
            
            # We should calculate the alignment score for both the 
            # sequences and the center residue using the HIJACK30 matrix
            score = nw.score_alignment("".join(seq1), "".join(seq2), 
                                      gap_open = 1000,
                                      gap_extend = 1000,
                                      matrix = matDir + 'HIJACK30')
            scoreCenter = nw.score_alignment(seq1Center, seq2Center,
                                            gap_open = 1000,
                                            gap_extend = 1000,
                                            matrix = matDir + 'HIJACK30')
            
            # The alignment scores should be scaled by 0.2
            scoreScaled = 0.2 * score
            scoreScaledCenter = 0.2 * scoreCenter
            
            # Correct the alignment score based on the center
            scoreCorrected = scoreScaled + scoreScaledCenter - 1000*0.2
            scores[i, j] = scoreCorrected
            
    return([scores, lengths])
    
