from pybicl.io import iterline, BedFile
from bx.intervals.cluster import ClusterTree
from pybicl.parser import Bed6Record
import sys, argparse

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="tmhmm_folder", required=True)
    base_group.add_argument("-r", "--ref", type=str, dest="reference", metavar="cds.bed12", required=True)
    base_group.add_argument("--out-genome", type=str, dest="out_genome", metavar="tmhmm.genome.bed12", required=True)
    base_group.add_argument("--out-cds", type=str, dest="out_cds", metavar="tmhmm.cds.bed6", required=True)
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
    f_out_genome = args.out_genome
    f_out_cds = args.out_cds
    f_ref = args.reference

    col_dict = {"outside": "200,50,50", "inside": "50,200,50", "TMhelix": "50,50,200"}
    cds_bed = BedFile(f_ref, "r")
    cds_dict = dict()
    for iso in cds_bed.load("isoform"):
        cds_dict[iso.name] = iso

    cds_li = list()
    with open(f_out_genome, "w") as f:
        for indx, line in enumerate(iterline(f_in)):
            if indx > 0:
                data = line.rstrip("\n").split()
                cds_name = data[0]
                pos = data[2]
                sidx = int(data[3]) - 1
                eidx = int(data[4])
                name = "{0}_{1}".format(cds_name, pos)
                cds_bed_rec = Bed6Record()
                cds_bed_rec.init_by_data(cds_name, sidx, eidx, pos, 0, ".")
                cds_li.append(cds_bed_rec)
                cds = cds_dict[cds_name]
                sub_cds = cds.slice(sidx * 3, eidx * 3).trans("bed12")
                sub_cds.score = 0
                sub_cds.itemRgb = col_dict[pos]
                sub_cds.name = name
                f.write(str(sub_cds) + "\n")
    cds_bed = BedFile(f_out_cds, "w")
    cds_bed.write(cds_li)

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
