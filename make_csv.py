import csv
import glob
import os

RECEPTOR_DIR = "/gs/hs0/tga-megadock/hayashi/immune/chain/"
LIGAND_DIR = "/gs/hs0/tga-megadock/hayashi/immune/chain/"
TARGET = "/gs/hs0/tga-megadock/hayashi/immune/result/"
DOCKING_OUTFILE = "/gs/hs0/tga-megadock/hayashi/data/result/immune.csv"

# Output: docking
# receptor name, ligand name, receptor size, ligand size, docking score...

def get_size(filename):
    size = 0
    with open(filename, "r") as f:
        for line in f:
            size = line[23:26].strip()
    return size

def get_score(filename):
    score_list = []
    with open(filename, "r") as f:
        for line in f:
            array = line.split("\t")
            if len(array) != 7:
                continue
            score = str(float(array[6].strip()))
            score_list.append(score)
    return score_list

def main():
    # get size of protein
    print("Start get protein size")
    # for receptor
    pdb_list = glob.glob(RECEPTOR_DIR + "*.pdb")
    protein_size_dct = {}
    for pdb_file in pdb_list:
        protein_name = os.path.basename(pdb_file).split(".")[0]
        #protein_name += "_r"
        protein_size = get_size(pdb_file)
        protein_size_dct[protein_name] = protein_size
    # for ligand
    """
    pdb_list = glob.glob(LIGAND_DIR + "*.pdb")
    for pdb_file in pdb_list:
        protein_name = os.path.basename(pdb_file).split(".")[0]
        #protein_name += "_l"
        protein_size = get_size(pdb_file)
        protein_size_dct[protein_name] = protein_size
    """
    print("Start out file parse")
    # parse out file
    result_list = []
    out_list = glob.glob(TARGET + "*.out")
    for outfile in out_list:
        print(outfile)
        file_name = os.path.basename(outfile)
        receptor = file_name.split("-")[0]
        #receptor += "_r"
        ligand = file_name.split("-")[1].split(".")[0]
        #ligand += "_l"
        score_list = get_score(outfile)
        score_list.insert(0, ligand)
        score_list.insert(0, receptor)
        score_list.insert(0, protein_size_dct[ligand])
        score_list.insert(0, protein_size_dct[receptor])
        result_list.append(score_list)
        
    print("Start write csv")
    # write csv
    with open(DOCKING_OUTFILE, "w") as f:
        writer = csv.writer(f)
        writer.writerows(result_list)

if __name__ == "__main__":
    main()
