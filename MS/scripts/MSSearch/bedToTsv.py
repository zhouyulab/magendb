#!/usr/bib/env python3

import os, sys, argparse

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.bed12", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.tsv", required=True)

    return(parser.parse_args(args))


def main():
    args = sys.argv[1:]

    args = parse_args(args)
    f_in = args.input
    f_out = args.output

    gene_db = {}
    with open(f_in, "r") as f:
        for line_data in f:
            line = line_data.strip("\n").split("\t")
            gene_id, pep_seq, position = line[3].split("|")
            PSMs = line[4]

            if gene_id not in gene_db:
                gene_db[gene_id] = {}
            gene_db[gene_id][pep_seq] = (position, PSMs)

    handle = open(f_out, "w")
    header = ["gene_id", "peptides", "positions", "PSMs"]
    handle.write("\t".join(header) + "\n")
    for gene_id in gene_db:
        pep_info = gene_db[gene_id]
        pep_all = [i for i in gene_db[gene_id]]
        position_all = [gene_db[gene_id][i][0].split(":")[1] for i in gene_db[gene_id]]
        psm_all = [gene_db[gene_id][i][1] for i in gene_db[gene_id]]
        
        info = [gene_id, ";".join(pep_all), ";".join(position_all), ";".join(psm_all)]
        handle.write("\t".join(info) + "\n")
    handle.close()




if __name__ == "__main__":
   main()  
