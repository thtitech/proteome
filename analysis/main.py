from protein import *
import sys
import csv
import glob
import numpy as np

DECOYNUM = 2000

def test_decoy():
    pdbfile = "/gs/hs0/tga-megadock/hayashi/immune/chain/1AU1_A.pdb"
    outfile = "/gs/hs0/tga-megadock/hayashi/immune/result/1A02_F-1AU1_A.out"
    tmpfile = "/gs/hs0/tga-megadock/hayashi/immune/test_decoy.txt"
    chain = parse_pdb(pdbfile)
    with open(outfile, "r") as f, open(tmpfile, "w") as ff:
        N, spacing = list(map(lambda x: float(x.strip()), next(f).split("\t")))
        rand_vec = list(map(lambda x: float(x.strip()), next(f).split("\t")[1:]))
        r_origin = list(map(lambda x: float(x.strip()), next(f).split("\t")[1:4]))
        l_origin = list(map(lambda x: float(x.strip()), next(f).split("\t")[1:4]))
        info = DockingInfo(N, spacing, rand_vec, r_origin, l_origin)
        # for test, top rank decoy is checked
        line = next(f).strip()
        array = list(map(lambda x: float(x.strip()), line.split("\t")))
        rotate_vec = array[0:3]
        parallel_vec = array[3:6]
        decoy = decoygen(chain, info, rotate_vec, parallel_vec)
        for res in decoy.res_list:
            #print(res.point)
            p = list(map(lambda x: str(x), res.point))
            ff.write(",".join(p) + "," + res.res_type + "," + str(res.res_num) + "\n")

def parse_info_file(info_file):
    chain_uni_dct = {}
    with open(info_file, "r") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            chain_uni_dct[row[3]] = row[2]
    return chain_uni_dct

def get_score(receptor, pred_recepotr, ligand, pred_ligand):
    # receptor and ligand are chain inastance
    # pred_receptor and pred_ligand are list of PredInterfaceResidue
    for r1 in receptor.res_list:
        for r2 in ligand.res_list:
            np.linalg.norm(np.array(r1.point) - np.array(r2.point))
    
    return 0

def main(pdb_dir, out_dir, seq_dir, info_file):

    log_f = open("/gs/hs0/tga-megadock/hayashi/proteome/analysis/progress_short.log", "w")
    log_f.write("Start\n")
    log_f.flush()


    chain_uni_dct = parse_info_file(info_file) # convert from chain to uniprot
    uni_pred_dct = {} # convert from uniprot to list of PredInterfaceResidue
    chain_list =[] # list of chain inastance
    
    pdb_file_list = glob.glob(pdb_dir + "*.pdb")
    seq_dir_list = glob.glob(seq_dir + "*")
    
    for seq_dir in seq_dir_list:
        uniprot_id = os.path.basename(seq_dir)
        uni_pred_dct[uniprot_id] = parse_pred_file(seq_dir + "/prediction.txt")

    for pdb_file in pdb_file_list:
        chain_id = os.path.basename(pdb_file).split(".")[0]
        pred_list = uni_pred_dct[chain_uni_dct[chain_id]]
        chain_list.append(parse_pdb(pdb_file, pred_list))

    for i, c1 in enumerate(chain_list):
        n = len(chain_list)
        log_f.write(c1.chain_id + ", " + str(i) + "/" + str(n) + "\n")
        log_f.flush()
        pred_c1 = uni_pred_dct[chain_uni_dct[c1.chain_id]]
        for c2 in chain_list:
            pred_c2 = uni_pred_dct[chain_uni_dct[c2.chain_id]]
            outfile = out_dir + c1.chain_id + "-" + c2.chain_id + ".out"
            with open(outfile, "r") as f:
                N, spacing = list(map(lambda x: float(x.strip()), next(f).split("\t")))
                rand_vec = list(map(lambda x: float(x.strip()), next(f).split("\t")[1:]))
                r_origin = list(map(lambda x: float(x.strip()), next(f).split("\t")[1:4]))
                l_origin = list(map(lambda x: float(x.strip()), next(f).split("\t")[1:4]))
                info = DockingInfo(N, spacing, rand_vec, r_origin, l_origin)
                # for test, top rank decoy is checked
                for i in range(DECOYNUM):
                    if i % 100 == 0:
                        print("Rank " + str(i))
                    line = next(f).strip()
                    array = list(map(lambda x: float(x.strip()), line.split("\t")))
                    rotate_vec = array[0:3]
                    parallel_vec = array[3:6]
                    decoy = decoygen(c2, info, rotate_vec, parallel_vec)
                    s = get_score(c1, pred_c1, decoy, pred_c2)
                    #print(c1.chain_id + "," + c2.chain_id + "," + str(i + 1) + "," + str(s))

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write("Invalid argment\n")
        sys.exit(0)
    tmp = list(map(lambda x: x + "/" if x[-1] != "/" else x, sys.argv[1:4]))
    pdb_dir = tmp[0]
    out_dir = tmp[1]
    seq_dir = tmp[2]
    info_file = sys.argv[4]
    main(pdb_dir, out_dir, seq_dir, info_file)
