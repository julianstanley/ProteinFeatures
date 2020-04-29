from src.utils import list_to_lookup
import copy

def get_protein_to_model(UniProtKB, protein_identifier, sequence, notes):
    '''
    Returns a formatted modeling form (standardized dictionary) 
    for the given protein sequence
    '''

    fasta_header = ">" + protein_identifier

    return({'UniProtKB' : UniProtKB,
        'fasta_header' : fasta_header,
        'fasta_sequence' : sequence,
        'notes' : notes})

def model_protein_segment(full_sequence, UniProtKB, segment_identifier, start, stop,
        notes = "protein fragment"):
    '''
    Create a formatted protein model for a segment of a larger protein

    full_sequence(str): fasta sequence of the full protein
    UniProtKB(str): UniProtKB identifier of the full protein
    segment_identifier(str): Some unique identifier for this segment
    start(int, [1, Inf]): The starting location of this segment
    end(start+1, Inf]): The ending location of this segment
    example: model_protein_segment("AKQQ", "P13807", "all_but_ala", 2, 4)
    '''

    # Format the identifier to ensure that
    # it does not have any spaces or hyphens
    segment_identifier = segment_identifier.replace(" ", "_")
    segment_identifier = segment_identifier.replace("-", "_")

    # Create the fasta header
    full_identifier = UniProtKB + "_" + segment_identifier \
    + "_" + str(start) + "_" + str(stop)

    # Take the start and stop segment
    fasta_sequence = full_sequence[start-1:stop]

    # Package it all up
    return(get_protein_to_model(UniProtKB, full_identifier, fasta_sequence, notes))

# Get pfams that are not covered by a crystal structure
def get_uncovered_identifiers(pfam_dicts, crystal_coverage):
    '''
    TODO: Implement
    '''

    uncovered_ids = []

    for pfam_dict in pfam_dicts:
        start = pfam_dict["Start"]
        stop = pfam_dict["Stop"]
        
        crystal_coverage_segment = crystal_coverage[start-1:stop]
        pfam_length = stop-start+1
        
        if sum(crystal_coverage_segment)/pfam_length < 0.95:
            pfam_accession = pfam_dict["Accession"]
            if pfam_accession not in uncovered_ids:
                uncovered_ids.append(pfam_accession)

    return(uncovered_ids)

def record_coverage(dictionary_list, total_length):
    '''
    dictionary_list must have "Start" and "Stop" keys
    corresponding to the start and stop locations of the entry 
    '''
    coverage = [0]*total_length

    for single_dict in dictionary_list:
        start = single_dict["Start"]
        stop = single_dict["Stop"]
        for loc in range(start-1, stop):
            coverage[loc] = 1

    return(coverage)

def should_model_flanking(pfam_dict, pfam_coverage):
    '''
    Should I model the flanking sequences of this pfam?
    '''

    # Get the starting and stopping positions of this pfam, and their 10bp up and 10bp down regions
    pfam_start = pfam_dict["Start"]
    pfam_stop = pfam_dict["Stop"]
    pfam_10down_start = pfam_start - 10
    pfam_10up_stop = pfam_stop + 10

    # Initalize whether we have room to expand up or downstream
    upstream = False
    downstream = False
    
    # Assuming we're not 10bp away from the start of the protein, 
    # we can expand downstream if there's no pfam coverage at that position    
    if pfam_10down_start >= 1:
        downstream = (pfam_coverage[pfam_10down_start-1] == 0)

    # Assuming we're not 10bp away from the end of the protein, 
    # we can expand downstream if there's no pfam coverage at that position   
    if pfam_10up_stop <= len(pfam_coverage): 
        upstream = (pfam_coverage[pfam_10up_stop-1] == 0)

    # We should model the flanking sequences if we have room to expand upstream OR downstream
    return(upstream or downstream)

def get_flanking_location(pfam, pfam_coverage):
    '''
    Gets the flanking regions of the given pfam
    Note: starting and ending locations start at 1, not 0
    '''

    pfam_start = pfam["Start"]
    pfam_stop = pfam["Stop"]

    # Get the starting position of the flanking region-----
    # Base case: if pfam starts at 1, then the flanking region also starts there
    if (pfam_start == 1):
        flanking_start = 0

    else:
        # Take all of the pfam coverage upstream of this pfam
        pfam_coverage_upstream = pfam_coverage[0:pfam_start-1]

        # If the upstream coverage has any pfams, then our flanking
        # region starts at the end of that pfam.
        if (1 in pfam_coverage_upstream):
            # Get the index of the last occurance of 1 in the list
            start_index = len(pfam_coverage_upstream) - 1 - pfam_coverage_upstream[::-1].index(1) 

            # Our start index should be the index after that last pfam
            flanking_start = start_index + 1

        # If there are no pfams upstream, we should start at the beginning
        else:
            flanking_start = 0

    # Get the ending position of the flanking region-----
    # Base case: if pfam ends at the end of the protein, then the flanking region goes
    # all the way to the end of the protein
    if (pfam_stop == len(pfam_coverage)):
        flanking_stop = len(pfam_coverage) - 1

    else:
        # Take all of the pfam coverage downstream of this pfam
        # Note: pfam_stop is indexed starting at 1
        pfam_coverage_downstream = pfam_coverage[pfam_stop:]

        # If the downstream coverage has any pfams, then our flanking
        # region starts at the beginning of that pfam
        if (1 in pfam_coverage_downstream):
            # Get the index of the first occurance of 1 in the list
            stop_index = pfam_coverage_downstream.index(1) + pfam_stop

            # Our stop index should be the index before that last pfam
            flanking_stop = stop_index - 1

        else:
            flanking_stop = len(pfam_coverage) - 1

    return((flanking_start + 1 , flanking_stop + 1))




# Get sequences 10bp up/downstream, combine pfams of same identifier, etc
def get_segments_of_interest(pfam_dicts):
    '''
    TODO: Implement
    NEED: Dictionary with UniProtKB --> length
    '''
    # Create a list of all new pfams to return
    pfam_dicts_new = []

    # Create a dictionary that goes {UniProtKB --> [list-of pfams]}
    lookup_uniprot_pfams = list_to_lookup(pfam_dicts, "UniProtKB")

    # Initalize a dictionary that will go {UniProtKB --> [pfam coverage]}
    lookup_uniprot_pfam_coverage = {}

    # Record pfam coverage for each protein.
    for UniProtKB in lookup_uniprot_pfams:
        sequence = lookup_uniprot_pfams[UniProtKB][0]["Sequence"]
        pfams = lookup_uniprot_pfams[UniProtKB]
        pfam_coverage = record_coverage(pfams, len(sequence))
        
        # Put the pfam coverage for this protein into the lookup table
        lookup_uniprot_pfam_coverage[UniProtKB] = pfam_coverage

    # Get segments to model from each pfam.
    for pfam_dict in pfam_dicts:
        # We want to model the original pfam, as long
        # as it isn't the entire protein (because we've already done that)
        start_loc = pfam_dict["Start"]
        stop_loc = pfam_dict["Stop"]
        sequence = pfam_dict["Sequence"]
        if not (start_loc == 1 and stop_loc == len(sequence)):
            pfam_dict["Notes"] = "Full pfam"
            pfam_dicts_new.append(copy.deepcopy(pfam_dict))

            # Get the pfam coverage of this protein
            uniprot_id = pfam_dict["UniProtKB"]
            pfam_coverage = lookup_uniprot_pfam_coverage[uniprot_id]

        # Do we need to model upstream and downstream?
        if(should_model_flanking(pfam_dict, pfam_coverage)):
            # Then we need to add the pfam + flanking regions to the queue
            start_loc, stop_loc = get_flanking_location(pfam_dict, pfam_coverage)

            # Don't model the entire protein (we already did that)
            if not (start_loc == 1 and stop_loc == len(pfam_dict["Sequence"])):
                new_pfam = pfam_dict
                new_pfam["Start"] = start_loc
                new_pfam["Stop"] = stop_loc
                new_pfam["Notes"] = "Pfam with flanking region"
                pfam_dicts_new.append(new_pfam)
        
    
    # We should also add big pfams when there are multiple pfams
    # with the same accession
    pfam_dicts_extra = []
    
    # Create a dictionary that goes {Accession --> [list-of pfams]}
    lookup_accession_pfams = list_to_lookup(pfam_dicts, "Accession")
    
    # Identify pfam accessions with duplicate entries, and concatenate them
    for accession in lookup_accession_pfams:
        pfams = lookup_accession_pfams[accession]
        if len(pfams) > 1:
            # Concatenate all pfams into one
            
            # initalize start and stop
            min_start = pfams[0]["Start"]
            max_stop = pfams[0]["Stop"]

            # Look through all pfams, picking min start and max stop
            for pfam in pfams:
                if pfam["Start"] < min_start:
                    min_start = pfam["Start"]

                if pfam["Stop"] > max_stop:
                    max_stop = pfam["Stop"]
            
            # Make a new pfam, add it to the list
            new_pfam = pfams[0]
            new_pfam["Start"] = min_start
            new_pfam["Stop"] = max_stop

            pfam_dicts_extra.append(new_pfam)


    # NOTE: repeating a loop here, but abstract this out later
        # Get segments to model from each pfam.
    for pfam_dict in pfam_dicts_extra:
        # We definitely want to model the original pfam
        pfam_dict["Notes"] = "Concatenated full pfam"
        pfam_dicts_new.append(pfam_dict)

        # Get the pfam coverage of this protein
        uniprot_id = pfam_dict["UniProtKB"]
        pfam_coverage = lookup_uniprot_pfam_coverage[uniprot_id]

        # Do we need to model upstream and downstream?
        if(should_model_flanking(pfam_dict, pfam_coverage)):
            # Then we need to add the pfam + flanking regions to the queue
            start_loc, stop_loc = get_flanking_location(pfam_dict, pfam_coverage)

            # Don't model the entire protein (we already did that)
            if not (start_loc == 1 and stop_loc == len(pfam_dict["Sequence"])):
                new_pfam = pfam_dict
                new_pfam["Start"] = start_loc
                new_pfam["Stop"] = stop_loc
                new_pfam["Notes"] = "Concatenated pfam with flanking region"

                pfam_dicts_new.append(new_pfam)

    return(pfam_dicts_new)


# We need to first make a list of all protein segments that we want to model
# Then we can model them
# What do we need to make that list? We need crystal structure and pfam
# information
# We need to know the pfam and the crystal information
def get_segments_to_model(UniProtKB, full_sequence, pfam_dicts, crystal_dicts):
    '''
    Get the segments that we should model for this protein
    '''
    # Initalize a list to keep track of what locations are covered by 
    # a crystal structure
    crystal_coverage = record_coverage(crystal_dicts, len(full_sequence))

    # Base case: if the crystal structures cover >95% of the structure,
    # then we don't need to model anything
    if sum(crystal_coverage)/len(full_sequence) > 0.95:
        return([])

    else:
        proteins_to_model = []
        
        # The whole protein should be modeled
        whole_protein_queue = get_protein_to_model(UniProtKB,
                UniProtKB + "_full_protein", full_sequence,
                "Full protein (incomplete crystal coverage)")
        proteins_to_model.append(whole_protein_queue)
        
        # Which pfam identifiers are incompletely covered?
        # Note: If one pfam from a group that shares an identifier 
        # is not covered, then all pfams in that group will be modeled
        uncovered_pfams = get_uncovered_identifiers(pfam_dicts, crystal_coverage)

        # Get the fully-enumerated list of pfams (including 10bp up/downstream
        # sequences, combined pfam families, etc.)
        # Effect: Dicts in pfam_dicts_full have an extra "Notes" entry
        # Describing 10bp up/down, combined, etc. 
        pfam_dicts_full = get_segments_of_interest(pfam_dicts)
        
        for pfam_dict in pfam_dicts_full:
            # Only worry about the pfam derived seq if its identifier
            # is in the uncovered_pfams list. If it's not, then we trust
            # that it has a crystal structure
            if pfam_dict["Accession"] in uncovered_pfams:
                # If this isn't covered, we need to add it to our list
                # of proteins to model
                UniProtKB = pfam_dict["UniProtKB"]
                segment_identifier = pfam_dict["Accession"] + "_" + pfam_dict["Identifier"]
                start = pfam_dict["Start"]
                stop = pfam_dict["Stop"]
                notes = pfam_dict["Notes"]

                # Call helper function to create a protein model with this segment
                segment_model = model_protein_segment(full_sequence, UniProtKB, 
                        segment_identifier, start, stop, notes)

                # Add this segment to the growing list of models
                proteins_to_model.append(segment_model)

        return(proteins_to_model)
