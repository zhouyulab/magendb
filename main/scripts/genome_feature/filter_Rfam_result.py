import argparse, sys
from pybicl.operation import IvalTree


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="Rfam.tbl", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="Rfam.filter.tbl", required=True)
    base_group.add_argument("--e-value", type=float, dest="e_value", metavar="e_value", default=1e-6)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = open(args.output, "w")
    e_value_cutoff = args.e_value

    ival_tree = IvalTree(True)

    with open(f_in, "r") as f:
        for indx, line in enumerate(f.readlines()):
            if indx in [0, 1]:
                f_out.write(line)
            else:
                if line.startswith("#"):
                    continue
                data = list(filter(lambda x: x, line.rstrip("\n").split(" ")))
                chrom = data[3]
                start = int(data[9])
                end = int(data[10])
                strand = data[11]
                score = float(data[16])
                evalue = float(data[17])
                if evalue > e_value_cutoff:
                    continue
                if strand == "+":
                    ival_tree.add((line, score), chrom, start, end, strand)
                else:
                    ival_tree.add((line, score), chrom, end, start, strand)

    clu_tree = ival_tree.trans2cluster()
    for ival in clu_tree.iter_all():
        overlap_recs = ival_tree.find(*ival)
        overlap_recs = sorted(overlap_recs, key=lambda x: x[1])
        line = overlap_recs[-1][0]
        f_out.write(line)
    f_out.close()


def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()