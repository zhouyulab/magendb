import argparse, sys
from pybicl.io import iterline
from collections import defaultdict
import numpy as np

def iter_data(file):
    for indx, line in enumerate(iterline(file)):
        if indx == 0:
            continue
        ori_prot1, ori_prot2, prot1, prot2, combined_score = line.rstrip().split("\t")
        yield prot1, prot2, float(combined_score)

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", nargs="*", type=str, dest="input", metavar="input.tsv", required=True)
    base_group.add_argument("--assemble", nargs="*", type=str, dest="assemble", metavar="assemble", required=True)
    base_group.add_argument("--uniq-locus-map", type=str, dest="uniq_locus", metavar="uniq_locus_map.tsv", required=True)
    base_group.add_argument("--out-prot", type=str, dest="out_prot", metavar="output.prot.tsv", required=True)
    base_group.add_argument("--out-locus", type=str, dest="out_locus", metavar="output.locus.tsv", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_uniq_locus = args.uniq_locus
    assemble_li = args.assemble
    f_out_prot = args.out_prot
    f_out_locus = args.out_locus

    assert len(f_in) == len(assemble_li)
    print("Loading score file ...")
    score_dict = defaultdict(list)
    for assemble, file in zip(assemble_li, f_in):
        for prot1, prot2, score in iter_data(file):
            prot1 = prot1.rstrip("_CDS")
            prot2 = prot2.rstrip("_CDS")
            score_dict[(assemble, prot1, prot2)].append(score)

    print("Building locus map ...")
    locus_map = dict()
    for indx, line in enumerate(iterline(f_uniq_locus)):
        if indx == 0:
            continue
        hash_id, assemble, cds_id = line.rstrip().split("\t")
        cds_id = cds_id.rstrip("_CDS")
        locus_map[(assemble, cds_id)] = hash_id
    locus_score_dict = defaultdict(list)

    print("Writing prot interactrion ...")
    with open(f_out_prot, "w") as f:
        header = "Assemble\tProtein1\tProtein2\tcombined_score\n"
        f.write(header)
        for (assemble, prot1, prot2), score_li in sorted(score_dict.items()):
            score = np.mean(score_li)
            data = [assemble, prot1, prot2, "{0:.1f}".format(score)]
            f.write("\t".join(data)+"\n")
            locus1 = locus_map[(assemble, prot1)]
            locus2 = locus_map[(assemble, prot2)]
            locus_score_dict[(locus1, locus2)].append(score)

    print("Writing locus interactrion ...")
    with open(f_out_locus, "w") as f:
        header = "Locus1\tLocus2\tcombined_score\n"
        f.write(header)
        for (locus1, locus2), score_li in sorted(locus_score_dict.items()):
            score = np.mean(score_li)
            data = [locus1, locus2, "{0:.1f}".format(score)]
            f.write("\t".join(data) + "\n")

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
