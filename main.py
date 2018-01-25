from chain import Chain
from kegg_converter import *
from uniprot_parser import *
from ppi_parser import *
import sys
import time
import itertools
from datetime import datetime

### Output Style
# pathway, kegg_id, uniprot_id, pdbid

### Output Style 2
# uniprotid1, uniprotid2

def main(args):
    # treat argments
    input_file = args[1]
    structure_output_file = args[2]
    interaction_output_file = args[3]
    # make kegg database
    kegg_converter = PathwayConverter()
    kegg_converter.make_database()
    # make interaction database
    # preppi_parser = PreppiParser()
    hint_parser = HintParser()
    #ppi_list = preppi_parser.get_ppi_list()
    #ppi_list.extend(hint_parser.get_ppi_list())
    ppi_list = hint_parser.get_ppi_list()
    
    with open(input_file, "r") as f, open(structure_output_file, "w") as sf, open(interaction_output_file, "w") as inf:
        sf.write(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + "\n")
        for line in f:
            pathway_id = line.strip()
            kegg_id_list = kegg_converter.get_gene_id_list(pathway_id)
            print("Start search pathway: " + pathway_id)
            uniprot_list_in_pathway = []
            count = 0
            all_chain_num = 0
            uniprot_num_with_no_chain = 0
            print("number of KEGG ID: " + str(len(kegg_id_list)))
            for kegg_id in kegg_id_list:
                count += 1
                if count % 10 == 0:
                    print(str(count) + "/" + str(len(kegg_id_list)))
                # get uniprot for each kegg id
                uniprot_id_list = kegg_to_uniprot(kegg_id)
                for uniprot_id in uniprot_id_list:
                    # uniprot_list_in_pathway.append(uniprot_id) if contain no chain ID
                    uniprot_converter = UniprotConverter(uniprot_id)
                    # use network connection
                    # convert to pdb id chain from
                    chain_list = uniprot_converter.get_chain_list()
                    if len(chain_list) == 0:
                        uniprot_num_with_no_chain += 1
                    else:
                        uniprot_list_in_pathway.append(uniprot_id)
                    all_chain_num += len(chain_list)
                    # if use filter, plaese change filter_chain()
                    chain_list = list(filter(lambda c: filter_chain(c), chain_list))
                    # filtering chain_list
                    chain_list = filter_chain_list(chain_list)
                    
                    for chain in chain_list:
                        column = []
                        column.append(pathway_id)
                        column.append(kegg_id)
                        column.append(uniprot_id)
                        column.append(chain.chain_name)
                        column.append(str(chain.start_res))
                        column.append(str(chain.end_res))
                        column.append(str(chain.resolution))
                        sf.write(",".join(column) + "\n")
                        
                time.sleep(1)
            print("number of uniprot id: " + str(len(uniprot_list_in_pathway)))
            print("number of uniprot id with no structure" + str(uniprot_num_with_no_chain))
            print("number of all chain id: " + str(all_chain_num))
            print("End search structure in pathway: " + pathway_id)
            print("Start search interaction in pathway: " + pathway_id)
            inf.write("#" + pathway_id + "\n")

            for a, b in list(itertools.combinations_with_replacement(set(uniprot_list_in_pathway),2)):
               if ((a, b) in ppi_list) or ((b, a) in ppi_list): 
                   inf.write(a + "," + b + "\n")
                        
            print("End search interaction in pathway: " + pathway_id)
            
def filter_chain(chain):
    # chain is Chain object, please refer chain.py
    # if you can filter resolution, resolution = -1 indicate this structure is not available with X-ray (ex. NMR or modeling ...)
    return ((chain.resolution != -1) and (chain.get_length() >= 30))

def filter_chain_list(chain_list):
    # sort by resolution -> get longer chain when overraping
    tmp_list = sorted(chain_list, key = lambda x: x.resolution)
    result = []
    for target in tmp_list:
        if len(result) == 0:
            result.append(target)
            continue
        if all(list(map(lambda x: not is_overrap_chain(target, x), result))):
            result.append(target)
    return result

def is_overrap_chain(c1, c2):
    # c1, c2 is chain instance
    # chain1: longer, chain2: shorter
    chain1 = c1
    chain2 = c2
    if chain1.get_length() < chain2.get_length():
        chain1 = c2
        chain2 = c1
    if (chain1.start_res > chain2.end_res) or (chain1.end_res < chain2.start_res):
        return False
    if (chain1.start_res <= chain2.start_res) and (chain1.end_res >= chain2.end_res):
        return True
    # set overrap
    overrap = 0
    if chain1.start_res >= chain2.start_res:
        overrap = chain2.end_res - chain1.start_res

    if chain1.end_res <= chain2.end_res:
        overrap = chain2.start_res - chain1.end_res

    if overrap / chain2.get_length() > 0.8:
        return True
    else:
        return False

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: python main.py [inputfile] [structure outputfile] [interaction outputfile]\n")
        sys.exit(0)
    main(sys.argv)
