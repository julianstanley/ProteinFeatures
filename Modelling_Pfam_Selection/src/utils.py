from itertools import compress

# A function to get an attribute value
def filter_on_attribute(element_list, attribute, value):
    # Get attributes of the typed passed in the attribute param
    # for each of the items in element_list
    attributes = [x.getAttribute(attribute) for x in element_list]
    
    # Filter to only get elements with the passed value
    filterParam = [x == value for x in attributes]
    return(list(compress(element_list, filterParam)))

# A function to turn a list of dictionaries into a lookup
# dictionary based on one of the common dictionary values
def list_to_lookup(dict_list, key):
    '''
    Takes a list of dictionaries and a common key
    and creates a lookup table such that 
    lookup_table[common key] --> dict
    '''
    
    lookup_table = {}

    for single_dict in dict_list:
        if single_dict[key] in lookup_table:
            lookup_table[single_dict[key]].append(single_dict)
        else:
            lookup_table[single_dict[key]] = [single_dict]

    return(lookup_table)
