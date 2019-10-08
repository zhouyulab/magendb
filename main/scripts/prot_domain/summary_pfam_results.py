from pybicl.io import iterline, BedFile
from bx.intervals.cluster import ClusterTree
import sys, argparse

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="pfam.out", required=True)
    base_group.add_argument("-r", "--ref", type=str, dest="reference", metavar="cds.bed12", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="pfam.bed12", required=True)
    return parser.parse_args(args)

def loc2region(li):
    clu = ClusterTree(0, 0)
    for x in li:
        clu.insert(x, x + 1, 0)
    for start, end, _ in clu.getregions():
        yield (start, end)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_ref = args.reference

    cds_bed = BedFile(f_ref, "r")
    cds_dict = dict()
    for iso in cds_bed.load("isoform"):
        cds_dict[iso.name] = iso
    with open(f_out, "w") as f:
        for line in iterline(f_in):
            if line:
                data = line.rstrip("\n").split()
                cds_name = data[0]
                sidx = int(data[1]) - 1
                eidx = int(data[2])
                name = "{0}_{1}:{2}({3})".format(cds_name, data[7], data[6], data[5])
                cds = cds_dict[cds_name]
                sub_cds = cds.slice(sidx * 3, eidx * 3).trans("bed12")
                sub_cds.score = 0
                sub_cds.name = name
                f.write(str(sub_cds) + "\n")

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
