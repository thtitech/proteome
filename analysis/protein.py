import os
import numpy as np
from numpy import sin, cos
import sys

class PredInterfaceResidue:
    def __init__(self, res_num, res_type, score):
        self.res_num = res_num
        self.res_type = res_type
        self.score = score

class Residue:
    # Residue is decsribed Ca atom
    def __init__(self, point, res_type, res_num):
        self.point = np.array(point)
        self.res_type = res_type
        self.res_num = res_num
    
class Chain:
    def __init__(self, chain_id):
        self.chain_id = chain_id
        self.res_list = []
        
    def add_res(self, res):
        self.res_list.append(res)

    def get_res(self, res_num):
        for res in self.res_list:
            if res.res_num == res_num:
                return res
        return None
        
class DockingInfo:
    def __init__(self, N, spacing, rand_vec, r_origin, l_origin):
        self.N = N
        self.spacing = spacing
        self.rand_vec = np.array(rand_vec)
        self.r_origin = np.array(r_origin)
        self.l_origin = np.array(l_origin)

def rotate(psi, theta, phi, point):
    r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi)
    r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi)
    r31 = sin(theta)*sin(phi)
    
    r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi)
    r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi)
    r32 = sin(theta)*cos(phi)
    
    r13 = sin(psi)*sin(theta)
    r23 = -cos(psi)*sin(theta)
    r33 = cos(theta)
    
    p = np.array(point)
    p = p.reshape(-1,)
    r = np.array([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])
    return np.dot(r, p).reshape(-1,)

def decoygen(chain, docking_info, rotate_vec, parallel_vec):
    decoy = Chain(chain.chain_id)
    N = docking_info.N
    spacing = docking_info.spacing
    for residue in chain.res_list:
        old_point = residue.point - docking_info.l_origin
        base_rotate_point = rotate(docking_info.rand_vec[0], docking_info.rand_vec[1], docking_info.rand_vec[2], old_point)
        lig_rotate_point = rotate(rotate_vec[0], rotate_vec[1], rotate_vec[2], base_rotate_point)
        t1 = parallel_vec[0] if (parallel_vec[0] < N/2) else parallel_vec[0] - N
        t2 = parallel_vec[1] if (parallel_vec[1] < N/2) else parallel_vec[1] - N
        t3 = parallel_vec[2] if (parallel_vec[2] < N/2) else parallel_vec[2] - N
        t = np.array([t1, t2, t3])
        new_point = lig_rotate_point - spacing * t + docking_info.r_origin
        new_res = Residue(new_point, residue.res_type, residue.res_num)
        decoy.add_res(new_res)
    return decoy

def parse_pdb(filename, pred_list):
    # file name must be ex. ~/1A12_A.pdb
    chain_id = os.path.basename(filename).split(".")[0]
    chain = Chain(chain_id)
    with open(filename, "r") as f:
        res_num = 0
        for line in f:
            if (line[0:4] != "ATOM") and (line[0:6] != "HETATM"):
                continue
            if line[13:15] == "CA":
                res_type = line[17:20]
                res_num += 1
                #res_num = int(line[22:26].strip())
                for pred in pred_list:
                    if (pred.res_num == res_num):
                        if (pred.res_type == res_type):
                            point = [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54])]
                            chain.add_res(Residue(point, res_type, res_num))
                        else:
                            print("Error in " + filename)
    return chain

def parse_pred_file(file_name):
    pred_list = []
    with open(file_name, "r") as f:
        for line in f:
            res_num = int(line[0:5].strip())
            res_type = line[8:11]
            score = float(line[16:].strip())
            pred_res = PredInterfaceResidue(res_num, res_type, score)
            pred_list.append(pred_res)
    return pred_list


