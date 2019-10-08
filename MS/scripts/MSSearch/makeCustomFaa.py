#!/usr/bin/env python3

import os, sys, argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type = str, dest = "input", metavar = "input.fasta", required = True)
    base_group.add_argument("-o", "--output", type = str, dest = "output", metavar = "output.fasta", required = True)
    return(parser.parse_args(args))


def main():
    args = sys.argv[1:]
    args = parse_args(args)
 
    f_input = args.input
    f_output = args.output

    with open(f_output, "w") as handle_f:
        with open(f_input, "r") as reader_f:
            for r in SeqIO.parse(reader_f, "fasta"):
                r_seq = r.seq.split("*")[0]
                if len(r_seq) >= 20:
                    r_custom = SeqRecord(r_seq, id = r.id, description = r.description)
                    SeqIO.write(r_custom, handle_f, "fasta")                
 


if __name__ == "__main__":
    main()
