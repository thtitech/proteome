import glob
import os

amino_code = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
    'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
    'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M',
    'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T',
    'TRP':'W', 'TYR':'Y', 'VAL':'V', 'HIS':'H',
    'ASX':'B', 'GLX':'Z', 'UNK':'K'
}

OUTFILE = "/gs/hs0/tga-megadock/hayashi/immune/align.txt"
PDBDIR = "/gs/hs0/tga-megadock/hayashi/immune/chain/"
SEQDIR = "/gs/hs0/tga-megadock/hayashi/immune/old/sequence/"
INFOFILE = "/gs/hs0/tga-megadock/hayashi/immune/immune_structure.txt"

pid_to_uni = {}
pid_list = []
with open(INFOFILE, "r") as f:
    next(f)
    for line in f:
        uni = line.strip().split(",")[2]
        pid = line.strip().split(",")[3]
        pid_to_uni[pid] = uni
        pid_list.append(pid)

uni_seq = {}
seq_dir = glob.glob(SEQDIR + "*")
for d in seq_dir:
    seq = ""
    with open(d + "/seq.txt", "r") as f:
        for line in f:
            if line[0] == ">":
                continue
            seq += line.strip()
        uni = os.path.basename(d)
        uni_seq[uni] = seq
    print(seq)

print("End make info")

with open(OUTFILE, "w") as f:
    N = len(pid_list)
    c = 0
    for pid in pid_list:
        c += 1
        print(str(c) + "/" + str(N) + "," + pid)
        seq = ""
        with open(PDBDIR + pid + ".pdb", "r") as ff:
            n = 0
            for line in ff:
                #if line[0:4] != "ATOM":
                    #continue
                tmp = int(line[22:26].strip())
                if tmp != n:
                    n = tmp
                    seq += amino_code[line[17:20]]
            
        base_seq = uni_seq[pid_to_uni[pid]]
        f.write(pid + "," + str(base_seq.find(seq)) + "\n")

