from html.parser import HTMLParser
from chain import Chain 

PREFIX = "http://www.uniprot.org/uniprot/"

class UniprotHTMLParser(HTMLParser):
    # UniprotのウェブページをパースしてPDB ID + chainを取得する
    def __init__(self, uniprot_id = ""):
        super().__init__()
        self.chain_name_list = []
        self.position_map = {}
        self.resolution_map = {}
        self.is_useful = False # 情報を吸い出したい場所か
        self.is_database_table = False # table内か
        self.is_finish = False
        self.count_tr = 0
        self.count_td = 0
        self.current_pdb_id = ""
        self.uniprot_id = uniprot_id

    def handle_starttag(self, tag, attrs):
        if self.is_finish:
            return
        if self.is_database_table and (tag == "tr"):
                self.count_tr += 1
        if self.is_database_table and (tag == "td"):
                self.count_td += 1
                if self.count_tr != 1:
                    self.is_useful = True
        for a in attrs:
            if a[0] == "class" and a[1] == "databaseTable STRUCTURE":
                # タグを見て、情報を吸い出す場所か判定
                self.is_database_table = True
                
    def handle_data(self, data):
        if self.is_finish:
            return
        if self.is_useful and (self.count_td == 1) and (len(data) == 4):
            self.current_pdb_id = data
            self.is_useful = False
        if self.is_useful and (self.count_td == 3):
            try:
                self.resolution_map[self.current_pdb_id] = float(data)
            except:
                return
        if self.is_useful and (self.count_td == 4):
            chains = data.split("/")
            for c in chains:
                self.chain_name_list.append(self.current_pdb_id + "_" + c)
        if self.is_useful and (self.count_td == 5):
            start, end = data.split("-")
            self.position_map[self.current_pdb_id] = (start, end)
        
    def handle_endtag(self, tag):
        if self.is_finish:
            return
        if self.is_database_table and (tag == "table"):
            self.is_finish = True
        elif self.is_database_table and (tag == "tr"):
            self.count_td = 0
    
    def get_url(self):
        return PREFIX + self.uniprot_id
    
    def get_chain_name_list(self):
        return self.chain_name_list
    
    def get_position_map(self):
        return self.position_map
    
    def get_resolution_map(self):
        return self.resolution_map

    def get_chain_list(self):
        result = []
        for chain_name in self.chain_name_list:
            pdb_id = chain_name.split("_")[0]
            resolution = self.resolution_map[pdb_id]
            start_res, end_res = self.position_map[pdb_id]
            result.append(Chain(pdb_id, resolution, start_res, end_res))
        return result
