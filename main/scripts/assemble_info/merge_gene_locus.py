import argparse, sys
from pybicl.operation import IvalTree
from pybicl.io import BedFile
from pybicl.utils import build_exon_overlap_isoform_cluster

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", nargs="*", type=str, dest="input", metavar="input.bed12", required=True)
    base_group.add_argument("-t", "--tag", nargs="*", type=str, dest="tag", metavar="tag", required=True)
    base_group.add_argument("-s", "--species", type=str, dest="species", metavar="species", required=True)
    base_group.add_argument("--gene-locus", type=str, dest="gene_locus", metavar="gene_locus.tsv", required=True)
    base_group.add_argument("--gene-locus-map", type=str, dest="gene_locus_map", metavar="gene_locus_map.tsv", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in_li = args.input
    tag_li = args.tag
    species = args.species
    f_gene_locus = args.gene_locus
    f_gene_locus_map = args.gene_locus_map

    assert len(f_in_li) == len(tag_li)
    bed_li = [BedFile(f, "r") for f in f_in_li]
    ival = IvalTree(True)
    for tag, bed in zip(tag_li, bed_li):
        for rec in bed.load("isoform"):
            rec.set_info(assemble = tag)
            ival.add_record(rec)
    merge_clu = ival.trans2cluster()

    f_locus = open(f_gene_locus, "w")
    f_locus_map = open(f_gene_locus_map, "w")
    locus_header = "GeneLocus\tChrom\tStart\tEnd\tStrand\tDomainAssemble\tDomainTid\n"
    locus_map_header = "GeneLocus\tAssemble\tTid\n"

    f_locus.write(locus_header)
    f_locus_map.write(locus_map_header)

    for loc in merge_clu.iter_all():
        chrom, start, end, strand = loc
        assert strand in ["+", "-"]
        overlap_rec = list(ival.find(*loc))
        for iso_clu in build_exon_overlap_isoform_cluster(overlap_rec):
            coding_iso_cds = list()
            for x in iso_clu:
                if x.thickStart != x.thickEnd:
                    cds = x.to_cds()
                    cds.set_info(assemble = x.get_info("assemble"))
            if coding_iso_cds:
                main_iso = sorted(coding_iso_cds, key=lambda x: x.exon_len())[-1]
            else:
                main_iso = sorted(iso_clu, key=lambda x: (x.blockCount, x.exon_len()))[-1]

            tmp_locus_start = min([x.chromStart for x in iso_clu])
            tmp_locus_end = max([x.chromEnd for x in iso_clu])

            if strand == "+":
                strand_str = "fwd"
            else:
                strand_str = "rev"
            locus_name = "{species}_{chrom}_{start}_{end}_{strand}".format(
                species=species, chrom=chrom, start=tmp_locus_start, end=tmp_locus_end, strand=strand_str
            )
            locus_text = "{locus}\t{chrom}\t{start}\t{end}\t{strand}\t{assemble}\t{tid}\n".format(
                locus=locus_name, chrom=chrom, start=tmp_locus_start, end=tmp_locus_end, strand=strand,
                assemble=main_iso.get_info("assemble"), tid=main_iso.name
            )
            f_locus.write(locus_text)
            for tmp_rec in iso_clu:
                locus_map_text = "{locus}\t{assemble}\t{tid}\n".format(
                    locus=locus_name, assemble=tmp_rec.get_info("assemble"), tid=tmp_rec.name
                )
                f_locus_map.write(locus_map_text)
    f_locus.close()
    f_locus_map.close()
    

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
