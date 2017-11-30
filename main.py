from chain import Chain
from kegg_converter import *
from uniprot_parser import *
from ppi_parser import *
import sys
import time

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
        for line in f:
            pathway_id = line.strip()
            kegg_id_list = kegg_converter.get_gene_id_list(pathway_id)
            print("Start search pathway: " + pathway_id)
            uniprot_list_in_pathway = []
            count = 0
            for kegg_id in kegg_id_list:
                count += 1
                if count % 10 == 0:
                    print(str(count) + "/" + str(len(kegg_id_list)))
                # get uniprot for each kegg id
                uniprot_id_list = kegg_to_uniprot(kegg_id)
                for uniprot_id in uniprot_id_list:
                    uniprot_list_in_pathway.append(uniprot_id)
                    uniprot_converter = UniprotConverter(uniprot_id)
                    # use network connection
                    # convert to pdb id chain from
                    chain_list = uniprot_converter.get_chain_list()
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
                        sf.write(",".join(column) + "\n")
                        
                time.sleep(1)
            print("End search structure in pathway: " + pathway_id)
            print("Start search interaction in pathway: " + pathway_id)
            inf.write("#" + pathway_id + "\n")
            for p1 in set(uniprot_list_in_pathway):
                for p2 in set(uniprot_list_in_pathway):
                    if ((p1, p2) in ppi_list) or ((p2, p1) in ppi_list):
                        inf.write(p1 + "," + p2 + "\n")
            print("End search interaction in pathway: " + pathway_id)
            print(uniprot_list_in_pathway)
            
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
