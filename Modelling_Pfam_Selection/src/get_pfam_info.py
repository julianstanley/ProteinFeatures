import requests
from xml.dom import minidom

# A function to get the uniprot information for a protein accession
def get_pfam_info(ID):
    # Initalize pfam data list
    pfam_data = []
    
    # Get pfam xml
    response = requests.get('http://pfam.xfam.org/protein/' + ID + '?output=xml')
    
    # Parse the xml into an object
    xmldoc = minidom.parseString(response.content)
    
    matches = xmldoc.getElementsByTagName('match')
    for match in matches:
        try:
            accession = str(match.getAttribute('accession'))
            identifier = str(match.getAttribute('id'))
            location = match.childNodes[1]
            start = int(location.getAttribute('start'))
            stop = int(location.getAttribute('end'))
            pfam_data.append({"UniProtKB" : ID,
                            "Accession" : accession,
                            "Identifier" : identifier,
                            "Start" : start,
                            "Stop" : stop})
        except Exception as e:
            print(e)
            continue
            
    return(pfam_data)       
