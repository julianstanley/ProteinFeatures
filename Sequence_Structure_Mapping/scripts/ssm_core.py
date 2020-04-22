#!/usr/bin/env python3
from glob import glob
from Bio import SeqIO
from warnings import warn
from os import path, makedirs
import pickle
from align import aligner
from datetime import date


def make_output_dirs(path_to_make):
    ''' Makes the given path, if it doesn't exist '''
    if not path.exists(path_to_make):
        makedirs(path_to_make)


def setup_output_dirs():
    ''' Sets up the output directories necessary for sequence structure mapping (ssm) '''
    dirs = ["./pickled_intermediates", "./output_maps/individual"]
    for path in dirs:
        make_output_dirs(path)


def build_sequence_dict(structure_fasta_path="./structure_fasta",
                        sequence_fasta_path="./sequence_fasta",
                        pickle_path="./pickled_intermediates"):
    ''' Builds a dictionary of structure files mapping to structure sequences and canonical sequences.
    All arguments have defaults. Finds structure files in structure_fasta_path

    Args:
        structure_fasta_path (str): Path to structure fasta files (corresponds to pdb files)
        sequence_fasta_path (str): Path to sequence fasta files (corresponds to canonical seqs)
        pickle_path (str): The path where we want to pickle this function's output
    Returns:
        (dict) that maps between structure file name and the corresponding structure and canonical sequences
    '''
    structure_to_sequences = {}
    for structure_fasta in glob(f'{structure_fasta_path}/*.fasta'):
        structure_fasta_basename = structure_fasta.split("/")[-1]
        uniprot_structure = structure_fasta_basename.replace(".fasta", "")
        up, structure = uniprot_structure.split("_", 1)
        structure_fasta_file = f'{structure_fasta_path}/{uniprot_structure}.fasta'
        sequence_fasta_file = f'{sequence_fasta_path}/{up}.fasta'

        if path.exists(structure_fasta_file) and path.exists(sequence_fasta_file):
            structure_fasta = str(
                [x for x in SeqIO.parse(structure_fasta_file, "fasta")][0].seq)
            sequence_fasta = str(
                [x for x in SeqIO.parse(sequence_fasta_file, "fasta")][0].seq)
            structure_to_sequences[uniprot_structure] = {"canonical": sequence_fasta,
                                                         "structure": structure_fasta}
        else:
            warn(
                f'Cannot find files: {structure_fasta_file} and {sequence_fasta_file}')

    with open(f"{pickle_path}/up_to_sequences.pickle", "wb") as f:
        pickle.dump(structure_to_sequences, f)

    return structure_to_sequences


def build_alignments(structure_to_sequences,
                     pickle_path="./pickled_intermediates"):
    ''' Given a map, align the canonical and structure-based sequences.
    Args:
        structure_to_sequences (dict, required. See build_sequence_dict). 
    '''
    alignments = {}
    for structure, sequences in structure_to_sequences.items():
        print(f"Aligning {structure}...")
        fasta_seq = sequences["canonical"].replace("U", "C")
        structure_seq = sequences['structure'].replace("U", "C")
        alignment = aligner(fasta_seq, structure_seq, method="glocal")[0]

        # If there are more than 5 mismatches, there may be a real gap.
        # Try a lower gap penalty
        if alignment.n_mismatches > 5:
            try:
                alignment = aligner(fasta_seq, structure_seq,
                                    method="global", gap_open=-5, gap_extend=-1)[0]
            except Exception as e:
                print(
                    "[NOTE] An exception is about to raise due to a failed global alignment.")
                print(f"[NOTE] {structure} had more than 5 mismatches with a semi-global",
                      f"alignment, but failed to generate any global alignments with gap_open = -5 and gap_extend = -1")
                raise

        alignments[structure] = alignment
        print(alignment)

    with open(f"{pickle_path}/alignments.pickle", "wb") as f:
        pickle.dump(alignments, f)

    return alignments


def get_structure_numbering(structure):
    ''' Gets the numbering from the given structure, identified by uniprot ID_structurename
    '''
    numbering_path = f"./structure_numbering/{structure}.pdb"
    if path.exists(numbering_path):
        with open(numbering_path) as f:
            return f.read().strip().split(",")
    else:
        warn(f'Numbering path {numbering_path} does not exist')
        return []


def make_alignment_map(alignment, numbering):
    # Goes from canonical sequence --> structure sequence
    alignment_map = {}

    # If there are no gaps, then everything matches up nicely.
    if alignment.n_gaps1 == 0 and alignment.n_gaps2 == 0:
        for i in range(0, alignment.end1 - alignment.start1):
            # Canonical sequence should start at 1 (i + 1), unless the first sequence alignment
            # doesn't start at 0, then we add accordingly.
            can_loc = i + 1 + alignment.start1

            # Structure sequence locations are given by the 'numbering' dict. So, just take the index
            # of each of these.
            struct_pos = i + alignment.start2
            struct_loc = int(numbering[struct_pos])

            # Throw into the dictionary
            alignment_map[can_loc] = struct_loc

    # If there are gaps, then we need to look at every character, so that we don't number gaps
    else:
        seq1 = alignment.seq1.decode("utf-8")
        seq2 = alignment.seq2.decode("utf-8")
        can_loc = alignment.start1 + 1
        struct_pos = alignment.start2
        for seq1_i, seq2_i in zip(seq1, seq2):
            if seq1_i == "-" or seq2_i == "-":
                alignment_map[can_loc] = "GAP"
            else:
                struct_loc = int(numbering[struct_pos])
                alignment_map[can_loc] = struct_loc

            if seq1_i != "-":
                can_loc += 1
            if seq2_i != "-":
                struct_pos += 1

    return alignment_map


def get_all_alignment_maps(alignments):
    maps = {}
    for structure, alignment in alignments.items():
        numbering = get_structure_numbering(structure)
        mapping = make_alignment_map(alignment, numbering)
        maps[structure] = {"map": mapping,
                           "canonical_seq": alignment.seq1.decode("utf-8"),
                           "structure_seq": alignment.seq2.decode("utf-8"),
                           "numbering": numbering}

    return maps


def make_alignment_csv(alignment_map, structure_to_sequences, output_dir="./output_maps"):
    with open(f"{output_dir}/map_{date.today()}.csv", "w") as full_outfh:
        header = "Sequence Name,Structure Name,Sequence Loc,Structure Loc\n"
        full_outfh.write(header)
        for structure, mapping in alignment_map.items():
            with open(f"{output_dir}/individual/{structure}.csv", "w") as individual_outfh:
                individual_outfh.write(header)
                canonical_seq = structure_to_sequences[structure]["canonical"].replace(
                    "U", "C")
                structure_seq = structure_to_sequences[structure]["structure"].replace(
                    "U", "C")
                numbering = mapping["numbering"]
                for canonical_loc, structure_loc in mapping["map"].items():
                    aa_canonical = canonical_seq[canonical_loc - 1]
                    if structure_loc == "GAP":
                        aa_structure = "GAP"
                    else:
                        aa_structure = structure_seq[numbering.index(
                            str(structure_loc))]
                    name_canonical, name_structure = structure.split("_", 1)
                    line = (f"{name_canonical}_{aa_canonical}{canonical_loc},",
                            f"{name_structure}.pdb_{aa_structure}{structure_loc},",
                            f"{canonical_loc},{structure_loc}\n")
                    full_outfh.write("".join(line))
                    individual_outfh.write("".join(line))


def main():
    setup_output_dirs()
    structure_to_seq = build_sequence_dict()
    alignments = build_alignments(structure_to_seq)
    maps = get_all_alignment_maps(alignments)
    make_alignment_csv(maps, structure_to_seq)


if __name__ == "__main__":
    main()
