#!/usr/bin/env python3

# Scripts influenced heavily by Michael J. Harms PDBTools

import sys, os

pdbtools.helper cmdline
pdbtools.data.common
pdbtools import seq

def parseCommandLine():
    """
    Parse command line.
    """

    global parser

    # ---------- Parse arguments --------------------

    options, args = parser.parse_args()

    # Generate list of pdb files on which to perform calculations, doing some
    # error checking along the way
    if len(args) < 1:
        err = "You must specify at least one pdb file!"
        parser.error(err)
    else:
        file_list, to_download = parseArgs(args)

    # Download missing pdb files
    if len(to_download) > 0:
        to_download.sort()
        if download.pdbDownload(to_download):
            file_list.extend(["%s.pdb" % f for f in to_download])
        else:
            err = "pdb could not be found on rcsb!"
            parser.error(err)

    # Remove duplicates from file_list by placing in dictionary
    file_dict = dict([(f,"") for f in file_list])
    file_list = list(file_dict.keys())
    file_list.sort()

    return file_list, options

def pdbSeq2Fasta(pdb,pdb_id="",chain="all",use_atoms=False):
    """
    Extract sequence from pdb file and write out in FASTA format.
    """

    # Grab sequences
    chain_dict, seq_type = pdbSeq(pdb,use_atoms)

    # Convert modified amino acids to their natural counterparts
    chain_dict = convertModifiedAA(chain_dict,pdb)

    # Determine which chains are being written out
    if chain == "all":
        chains_to_write = list(chain_dict.keys())
        chains_to_write.sort()
    else:
        if chain in list(chain_dict.keys()):
            chains_to_write = [chain]
        else:
            err = "Chain \"%s\" not in pdb!" % chain
            raise PdbSeqError(err)

    # Convert sequences to 1-letter format and join strings
    for c in chains_to_write:
        for aa_index, aa in enumerate(chain_dict[c]):
            try:
                chain_dict[c][aa_index] = AA3_TO_AA1[aa]
            except KeyError:
                chain_dict[c][aa_index] = "X"

    out = []
    for c in chains_to_write:
        out.append(">%s%s_%s" % (pdb_id,c,seq_type))

        # Write output in lines 80 characters long
        seq_length = len(chain_dict[c])
        num_lines = seq_length // 80

        for i in range(num_lines+1):
            out.append("".join([aa for aa in chain_dict[c][80*i:80*(i+1)]]))
        out.append("".join([aa for aa in chain_dict[c][80*(i+1):]]))
        if out[-1] == "":
            out.pop(-1)


    return "\n".join(out)

def pdbSeq(pdb,use_atoms=False):
    """
    Parse the SEQRES entries in a pdb file.  If this fails, use the ATOM
    entries.  Return dictionary of sequences keyed to chain and type of
    sequence used.
    """

    # Try using SEQRES
    seq = [l for l in pdb if l[0:6] == "SEQRES"]
    if len(seq) != 0 and not use_atoms:
        seq_type = "SEQRES"
        chain_dict = dict([(l[11],[]) for l in seq])
        for c in list(chain_dict.keys()):
            chain_seq = [l[19:70].split() for l in seq if l[11] == c]
        for x in chain_seq:
                chain_dict[c].extend(x)

    # Otherwise, use ATOM
    else:

        seq_type = "ATOM  "

        # Check to see if there are multiple models.  If there are, only look
        # at the first model.
        models = [i for i, l in enumerate(pdb) if l.startswith("MODEL")]
        if len(models) > 1:
            pdb = pdb[models[0]:models[1]]

        # Grab all CA from ATOM entries, as well as MSE from HETATM
        atoms = []
        for l in pdb:
            if l[0:6] == "ATOM  " and l[13:16] == "CA ":

                # Check to see if this is a second conformation of the previous
                # atom
                if len(atoms) != 0:
                    if atoms[-1][17:26] == l[17:26]:
                        continue

                atoms.append(l)
            elif l[0:6] == "HETATM" and l[13:16] == "CA " and l[17:20] == "MSE":

                # Check to see if this is a second conformation of the previous
                # atom
                if len(atoms) != 0:
                    if atoms[-1][17:26] == l[17:26]:
                        continue

                atoms.append(l)

        chain_dict = dict([(l[21],[]) for l in atoms])
        for c in list(chain_dict.keys()):
            chain_dict[c] = [l[17:20] for l in atoms if l[21] == c]

    return chain_dict, seq_type

def main():
    """
    Function to execute if called from command line.
    """
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                      long_flag="chain",
                      action="store",
                      default="all",
                      help="chain to select",
                      nargs=1,
                      type=str)
    cmdline.addOption(short_flag="a",
                      long_flag="atomseq",
                      action="store_true",
                      default=False,
                      help="use ATOM sequence, not SEQRES")


    file_list, options = cmdline.parseCommandLine()

    # Extract sequence data
    for pdb_file in file_list:

        pdb_id = os.path.split(pdb_file)[-1][:-4]

        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        sequence = pdbSeq2Fasta(pdb,pdb_id,options.chain,options.atomseq)

        print(sequence)



if __name__ == "__main__":
    main()
