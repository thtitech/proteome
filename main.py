from chain import Chain
from kegg_converter import *
from uniprot_parser import *
import sys
import time

### Output Style
# pathway, kegg_id, uniprot_id, pdbid

def main(args):
    input_file = args[1]
    output_file = args[2]
    kegg_converter = PathwayConverter()
    kegg_converter.make_database()
    with open(input_file, "r") as f, open(output_file, "w") as ff:
        for line in f:
            pathway_id = line.strip()
            kegg_id_list = kegg_converter.get_gene_id_list(pathway_id)
            print("Start search pathway: " + pathway_id)
            for kegg_id in kegg_id_list:
                uniprot_id_list = kegg_to_uniprot(kegg_id)
                for uniprot_id in uniprot_id_list:
                    uniprot_converter = UniprotConverter(uniprot_id)
                    # use network connection
                    chain_list = uniprot_converter.get_chain_list()
                    # if use filter, plaese change filter_chain()
                    chain_list = list(filter(lambda c: filter_chain(c), chain_list))
                    for chain in chain_list:
                        column = []
                        column.append(pathway_id)
                        column.append(kegg_id)
                        column.append(uniprot_id)
                        column.append(chain.chain_name)
                        ff.write(",".join(column) + "\n")

def filter_chain(chain):
    # chain is Chain object, please refer chain.py
    # if you can filter resolusion, resolution = -1 indicate this structure is not available with X-ray (ex. NMR or modeling ...)
    return True

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python main.py [inputfile] [outputfile]\n")
        sys.exit(0)
    main(sys.argv)



    
