import argparse, sys
from pybicl.io import BedFile, iterline


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.out", required=True)
    base_group.add_argument("-r", "--ref", type=str, dest="ref", metavar="reference.bed12", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.bed12", required=True)
    base_group.add_argument("--min-loop-len", type=int, dest="min_loop_len", metavar="min_loop_len", default=1)
    base_group.add_argument("--min-rg4-len", type=int, dest="min_rg4_len", metavar="min_rg4_len", default=3)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_ref = args.ref
    f_output = args.output
    min_loop_len = args.min_loop_len
    min_rg4_len = args.min_rg4_len

    utr_bed = BedFile(f_ref, "r")
    utr_recs = dict()
    for rec in utr_bed.load("isoform"):
        utr_recs[rec.name] = rec
    res_rec = list()
    idx = 0
    for line in iterline(f_in):
        name, rg4_indx, sidx1, sidx2, sidx3, sidx4, rg4_len, score, seq = line.rstrip("\n").split("\t")
        sidx1 = int(sidx1)
        sidx2 = int(sidx2)
        sidx3 = int(sidx3)
        sidx4 = int(sidx4)
        rg4_len = int(rg4_len)
        if rg4_len < min_rg4_len:
            continue
        if (sidx2 - sidx1 - rg4_len - min_loop_len) < 0:
            continue
        if (sidx3 - sidx2 - rg4_len - min_loop_len) < 0:
            continue
        if (sidx4 - sidx3 - rg4_len - min_loop_len) < 0:
            continue
        idx += 1
        rec = utr_recs[name]
        rg4_region = rec.slice(sidx1, sidx4 + rg4_len, name="{0}_rg4_{1}".format(name, idx), strand=True)
        res_rec.append(rg4_region.trans("bed12"))
    out_bed = BedFile(f_output, "w")
    out_bed.write(res_rec)

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
