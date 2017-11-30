import urllib
import requests
import time

DBURL = "http://rest.kegg.jp/link/pathway/hsa"
CONVERT_URL = "http://www.uniprot.org/uploadlists/"

class PathwayConverter:
    # Pathway に含まれるタンパク質を取得する
    def __init__(self):
        self.database = {}

    def make_database(self):
        r = requests.get(DBURL)
        for line in r.text.split("\n"):
            array = line.split("\t")
            if len(array) != 2:
                continue
            gene_id, pathway_id = array
            pathway_id = pathway_id.split(":")[1]
            if pathway_id not in self.database:
                self.database[pathway_id] = []
            self.database[pathway_id].append(gene_id)

    def get_gene_id_list(self, pathway_id):
        return self.database.get(pathway_id, [])
            
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
        
def kegg_list_to_uniprot(kegg_id_list):
    result = {}
    for kegg_id in kegg_id_list:
        res = kegg_to_uniprot(kegg_id)
        result[kegg_id] = res
        time.sleep(0.1)
    return result

