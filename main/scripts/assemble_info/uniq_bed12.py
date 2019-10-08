import argparse, sys
from pybicl.io import BedFile

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", nargs="*", type=str, dest="input", metavar="input.bed12", required=True)
    base_group.add_argument("-t", "--tag", nargs="*", type=str, dest="tag", metavar="tag", required=True)
    base_group.add_argument("--uniq", type=str, dest="uniq_rec", metavar="uniq.bed12", required=True)
    base_group.add_argument("--uniq-map", type=str, dest="uniq_map", metavar="uniq_map.tsv", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in_li = args.input
    tag_li = args.tag
    f_uniq_rec = args.uniq_rec
    f_uniq_map = args.uniq_map

    key_li = set()
    uniq_rec_li = list()
    id_set = set()
    with open(f_uniq_map, "w") as f:
        header = "HashID\tAssemble\tName\n"
        f.write(header)
        for tag, f_in in zip(tag_li, f_in_li):
            bed = BedFile(f_in, "r")
            for rec in bed.load():
                key = (rec.chrom, rec.chromStart, rec.chromEnd, rec.strand, rec.thickStart, rec.thickEnd, tuple(rec.iterblock()))
                tmp_id = "{0}_{1}".format(rec.chrom, abs(hash(key)))
                data = "{0}\t{1}\t{2}\n".format(tmp_id, tag, rec.name)
                f.write(data)
                if key not in key_li:
                    key_li.add(key)
                    assert tmp_id not in id_set
                    id_set.add(tmp_id)
                    rec.name = tmp_id
                    uniq_rec_li.append(rec)
    out_bed = BedFile(f_uniq_rec, "w")
    out_bed.write(uniq_rec_li)

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
