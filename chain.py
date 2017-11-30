class Chain:
    def __init__(self, chain_name, resolution, start_res = 0, end_res = 0):
        self.chain_name = chain_name
        self.resolution = resolution
        self.start_res = start_res
        self.end_res = end_res

    def get_pdb_id(self):
        return self.chain_name.split("_")[0]

    def get_length(self):
        return (self.end_res - self.start_res + 1)

    def __repr__(self):
        return self.chain_name
