import argparse, sys
from pybicl.io import iterline
from collections import defaultdict

def load_ref_name(f):
    name_li = list()
    for line in iterline(f):
        data = line.rstrip().split("\t")
        name_li.append(data[1])
    return name_li

def iter_coliner_pair(f, string_li, plant_li):
    for line in iterline(f):
        data = line.rstrip().split("\t")
        name1 = data[1]
        name2 = data[2]
        if name1 in string_li:
            if name2 in plant_li:
                yield name1, name2
        else:
            assert name1 in plant_li
            if name2 in string_li:
                yield name2, name1

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--coliner", type=str, dest="coliner", metavar="coliner.tsv", required=True)
    base_group.add_argument("--ppi", type=str, dest="ppi", metavar="ppi.tsv", required=True)
    base_group.add_argument("--string-ref", type=str, dest="string_ref", metavar="string_ref.tsv", required=True)
    base_group.add_argument("--plant-ref", type=str, dest="plant_ref", metavar="plant_ref.tsv", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.tsv", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_coliner = args.coliner
    f_ppi = args.ppi
    f_string_ref = args.string_ref
    f_plant_ref = args.plant_ref
    f_out = args.output

    string_name_li = load_ref_name(f_string_ref)
    plant_name_li = load_ref_name(f_plant_ref)
    trans_map = defaultdict(set)
    for string_name, plant_name in iter_coliner_pair(f_coliner, string_name_li, plant_name_li):
        trans_map[string_name].add(plant_name)

    with open(f_out, "w") as f:
        header = "ori_prot1\tori_prot2\tprot1\tprot2\tcombined_score\n"
        f.write(header)
        for indx, line in enumerate(iterline(f_ppi)):
            if indx == 0:
                continue
            ori_prot1, ori_prot2, score = line.rstrip().split(" ")
            for prot1 in sorted(trans_map[ori_prot1]):
                for prot2 in sorted(trans_map[ori_prot2]):
                    if prot1 == prot2:
                        continue
                    elif prot1 < prot2:
                        data = [ori_prot1, ori_prot2, prot1, prot2, score]
                    else:
                        data = [ori_prot2, ori_prot1, prot2, prot1, score]
                    f.write("\t".join(data)+"\n")

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()