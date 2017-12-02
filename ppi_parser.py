import csv

PREPPI1 = "preppi_1.txt"
PREPPI2 = "preppi_2.txt"
HINT = "hint.txt"

class PreppiParser:
    def get_ppi_list(self):
        result = []
        result.append(self.parse_file(PREPPI1))
        result.append(self.parse_file(PREPPI2))
        return result
    
    def parse_file(self, filename):
        result = []
        with open(filename, "r") as f:
            for line in f:
                if line[0] == "#":
                    continue
                db_id, uniprot_pair, dbs, pubmeds = line.strip().split("\t")
                p1, p2 = uniprot_pair.split("|")
                result.append((p1, p2))
        return result

class HintParser:
    def get_ppi_list(self):
        result = []
        with open(HINT, "r") as f:
            reader = csv.reader(f, delimiter = "\t")
            next(reader)
            for row in reader:
                result.append((row[0], row[1]))
        return result
