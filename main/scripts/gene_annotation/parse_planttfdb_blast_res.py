from pybicl.io import BlastXML
import argparse
import sys
import re

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", nargs="*", type=str, dest="input", metavar="blast.out", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="annotation.tsv", required=True)
    base_group.add_argument("--e-val", type=float, dest="e_val", metavar="e_value", default=1e-5)
    base_group.add_argument("--hsp-ratio", type=float, dest="hsp_ratio", metavar="hsp_ratio", default=0)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    lim_eval = args.e_val
    lim_hsp_ratio = args.hsp_ratio
    tf_family_pattern = re.compile("\|([\w\-/]+)\|")
    hist = list()
    with open(f_out, "w") as f:
        f.write("TfFamily\tName\tEvalue\tHspScore\tHspRatio\n")
        for blastout in f_in:
            blast_file = BlastXML(blastout, "r")
            blast_root = blast_file.load()
            for query in blast_root.iter_query():
                query_name = query.name
                hit_iter = query.iter_hit()
                while True:
                    try:
                        hit = next(hit_iter)
                    except StopIteration:
                        break
                    hsp = hit.hsp_list[0]
                    hit_id = hit.id
                    e_val = hsp.eval
                    score = hsp.score
                    if lim_eval < e_val:
                        continue
                    query_len = query.len
                    hit_len = hit.len
                    hsp_len = hsp.align_len
                    hsp_ratio = min(1, hsp_len / min(query_len, hit_len))
                    if hsp_ratio < lim_hsp_ratio:
                        continue
                    tf_family = re.findall(tf_family_pattern, hit_id)[0]
                    key = (tf_family, query_name)
                    if key in hist:
                        continue
                    hist.append(key)
                    f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                        tf_family, query_name, e_val, score, hsp_ratio
                    ))

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
