import requests
from xml.dom import minidom
from src.utils import filter_on_attribute
import requests

# A function to get the uniprot information for a protein accession
def get_uniprot_info(ID):
    # Get uniprot xml
    response = requests.get('https://www.uniprot.org/uniprot/' + ID + '.xml')
    
    # Parse the xml into an object
    xmldoc = minidom.parseString(response.content)
    
    # Get all PDB-type dbReference tags
    PDB_refs = filter_on_attribute(
    element_list = xmldoc.getElementsByTagName('dbReference'),
    attribute = "type",
    value = "PDB")
    
    # Add all X-ray data to store in a dictionary
    PDB_data = []
    for structure in PDB_refs:
        model_name = str(structure.getAttribute('id'))
        nodes = structure.childNodes
        try:
            method = nodes[1].getAttribute('value')
            some_crystalization = False
            if method == "X-ray":  
                resolution = float(nodes[3].getAttribute('value'))
                start_stop = nodes[5].getAttribute('value').split("=")[1].split("-")   
                start = int(start_stop[0])
                stop = int(start_stop[1])

                PDB_data.append({'UniProtKB' : ID, 
                                'Model Name' : model_name,
                                'Resolution' : resolution,
                                'Start' : start,
                                'Stop' : stop})
        except:
            continue

    # Get the protein sequence while we're at it
    sequence_nodes = xmldoc.getElementsByTagName('sequence')
    sequence_node = sequence_nodes[len(sequence_nodes)-1].childNodes[0]
    sequence = str(sequence_node.data).strip().replace("\n", "")

    return(PDB_data, sequence)       
