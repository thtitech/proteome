import requests
from xml.etree.ElementTree import fromstring

pdb_id = '4hhb.A'
pdb_mapping_url = 'http://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment'
uniprot_url = 'http://www.uniprot.org/uniprot/{}.xml'

def get_uniprot_accession_id(response_xml):
    root = fromstring(response_xml)
    return next(
        el for el in root.getchildren()[0].getchildren()
        if el.attrib['dbSource'] == 'UniProt'
    ).attrib['dbAccessionId']

def get_uniprot_protein_name(uniport_id):
    uinprot_response = requests.get(
        uniprot_url.format(uniport_id)
    ).text
    return fromstring(uinprot_response).find(
        './/{http://uniprot.org/uniprot}recommendedName/{http://uniprot.org/uniprot}fullName'
    ).text

def map_pdb_to_uniprot(pdb_id):
    pdb_mapping_response = requests.get(
        pdb_mapping_url, params={'query': pdb_id}
    ).text
    uniprot_id = get_uniprot_accession_id(pdb_mapping_response)
    uniprot_name = get_uniprot_protein_name(uniprot_id)
    return {
        'pdb_id': pdb_id,
        'uniprot_id': uniprot_id,
        'uniprot_name': uniprot_name
    }

print (map_pdb_to_uniprot(pdb_id))
