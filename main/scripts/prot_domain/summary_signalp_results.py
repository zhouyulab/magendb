from pybicl.io import iterline, BedFile
from pybicl.parser import Bed6Record
import sys, argparse, os, re

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="signalp.out", required=True)
    base_group.add_argument("-r", "--ref", type=str, dest="reference", metavar="cds.bed12", required=True)
    base_group.add_argument("--out-genome", type=str, dest="out_genome", metavar="signalp.genome.bed12", required=True)
    base_group.add_argument("--out-cds", type=str, dest="out_cds", metavar="signalp.cds.bed6", required=True)
    return parser.parse_args(args)

def iter_signalp_rec(f):
    re_compile_bool = re.compile("Name=\S+\s+SP='(\w+)'\s+")
    re_compile_data = re.compile("Name=\S+\s+SP='\w+'\s+Cleavage site between pos. (\d+) and \d+:\s+(\S+)")
    cds_name = None
    pos = None
    site = None
    for line in iterline(f):
        if line[0] == ">":
            cds_name = line.rstrip().lstrip(">")
            continue
        data = re.findall(re_compile_bool, line)
        if data and data[0] == "YES":
            pos, site = re.findall(re_compile_data, line)[0]
            pos = int(pos) - 1
            yield cds_name, pos, site
        

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out_genome = args.out_genome
    f_out_cds = args.out_cds
    f_ref = args.reference

    cds_bed = BedFile(f_ref, "r")
    cds_dict = dict()
    for iso in cds_bed.load("isoform"):
        cds_dict[iso.name] = iso

    cds_li = list()
    with open(f_out_genome, "w") as f:
        for cds_name, pos, site in iter_signalp_rec(f_in):
            cds = cds_dict[cds_name]
            cds_bed_rec = Bed6Record()
            cds_bed_rec.init_by_data(cds_name, 0, pos, site, 0, ".")
            cds_li.append(cds_bed_rec)
            sub_cds = cds.slice(0, pos * 3).trans("bed12")
            sub_cds.score = 0
            sub_cds.name = "{0}:{1}".format(cds_name, site)
            f.write(str(sub_cds) + "\n")
    cds_bed = BedFile(f_out_cds, "w")
    cds_bed.write(cds_li)

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
