import urllib
import requests

DBURL = "http://rest.kegg.jp/link/pathway/hsa"
CONVERT_URL = "http://www.uniprot.org/uploadlists/"

class PathwayConverter:
    # Pathway に含まれるタンパク質を取得する
    def __init__(self):
        self.database = {}

    def make_database(self):
        r = requests.get(DBURL)
        for line in r.text.split("\n"):
            gene_id, pathway_id = line.split("\t")
            pathway_id = pathway_id.split(":")
            if pathway_id not in self.database:
                self.database[pathway_id] = []
            self.database[pathway_id].append(gene_id)

    def get_gene_id_list(self, pathway_id):
        if pathway_id in self.database:
            return self.database[pathway_id]
        else:
            return []
            
def kegg_to_uniprot(kegg_id):
    # KEGG のhsa:~のidをuniprotに変換してくれる
    result = []
    params = {"from": "KEGG_ID", "to": "ACC", "format": "tab"}
    params["query"] = kegg_id
    r = requests.get(CONVERT_URL, params)
    for line in r.text.split("\n"):
        a = line.split("\t")
        if a[0].lower() == kegg_id:
            result.append(a[2])
    return result
        
