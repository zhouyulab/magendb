from pybicl.io import BlastXML
import argparse
import sys

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", nargs="*", type=str, dest="input", metavar="blast.out", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="annotation.tsv", required=True)
    base_group.add_argument("--e-val", type=float, dest="e_val", metavar="e_value", default=1e-5)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    lim_eval = args.e_val
    text = "{qname}\t{hit_id}\t{hit_acc}\t{identity}\t{align_len}\t{q_start}\t{q_end}\t{t_start}\t{t_end}\t{e_val}\t{score}\t{hit_def}\n"
    with open(f_out, "w") as f:
        f.write("#QueryName\tHitId\tHitAccession\tIdentity\tAlignLength\tQueryStart\tQueryEnd\tHitStart\tHitEnd\tE_Value\tScore\tFunction\n")
        for blastout in f_in:
            blast_file = BlastXML(blastout, "r")
            blast_root = blast_file.load()
            for query in blast_root.iter_query():
                query_name = query.name
                for hit in query.iter_hit():
                    hsp = hit.hsp_list[0]
                    hit_id = hit.id
                    hit_accession = hit.accession
                    identity = hsp.identity
                    align_len = hsp.align_len
                    q_start = hsp.query_from
                    q_end = hsp.query_to
                    t_start = hsp.hit_from
                    t_end = hsp.hit_to
                    e_val = hsp.eval
                    score = hsp.score
                    hit_def = hit.func
                    if lim_eval < e_val:
                        continue
                    f.write(text.format(
                        qname=query_name, hit_id=hit_id, hit_acc=hit_accession, identity=identity,
                        align_len=align_len, q_start=q_start, q_end=q_end, t_start=t_start,
                        t_end=t_end, e_val=e_val, score=score, hit_def=hit_def
                    ))

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
