from pdb_to_uniprot_mapper import *
import sys

pid_list = set()
with open("/gs/hs0/tga-megadock/hayashi/benchmark2/protein_info.txt", "r") as f:
    for line in f:
        array = line.split(",")
        pid_list.add(array[1].strip().replace("_", "."))
        pid_list.add(array[2].strip().replace("_", "."))

with open("/gs/hs0/tga-megadock/hayashi/benchmark2/seq_list.txt", "w") as f:
    N = len(pid_list)
    c = 1
    for pid in pid_list:
        print(str(c) + "/" + str(N))
        c+=1
        try:
            d = map_pdb_to_uniprot(pid)
        except Exception as e:
            f.write(pid.replace(".", "_") + "," + "\n")
        else:
            f.write(pid.replace(".", "_") + "," + d["uniprot_id"] + "\n")
    
