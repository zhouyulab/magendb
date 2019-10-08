import argparse, sys
from pybicl.io import BedFile, iterline
from pybicl.operation import IvalTree
import re
from collections import defaultdict

class GenomeInfo(object):
    def __init__(self, assemble, f_cds_ref, f_uniq_cds_map, f_uniq_locus_map):
        self.assemble = assemble
        self.cds_bed = BedFile(f_cds_ref, "r")
        self.cds_ival = IvalTree(True)
        self.cds_dict = dict()
        self.uniq_cds_map = dict()
        self.uniq_locus_map = dict()
        self.load_cds()
        self.load_uniq_cds_map(f_uniq_cds_map)
        self.load_uniq_locus_map(f_uniq_locus_map)

    def load_cds(self):
        for rec in self.cds_bed.load():
            self.cds_ival.add_record(rec)
            self.cds_dict[rec.name.rstrip("_CDS")] = rec

    def load_uniq_cds_map(self, f_uniq_cds_map):
        for indx, line in enumerate(iterline(f_uniq_cds_map)):
            if indx == 0:
                continue
            hash_id, assemble, cds_id = line.rstrip().split("\t")
            if assemble != self.assemble:
                continue
            cds_id = cds_id.rstrip("_CDS")
            self.uniq_cds_map[cds_id] = hash_id

    def load_uniq_locus_map(self, f_uniq_locus_map):
        for indx, line in enumerate(iterline(f_uniq_locus_map)):
            if indx == 0:
                continue
            locus_id, assemble, tid = line.rstrip().split("\t")
            if assemble != self.assemble:
                continue
            self.uniq_locus_map[tid] = locus_id


    def is_in_genome(self, name):
        return name in self.cds_dict.keys()

    def get_cds(self, name):
        return self.cds_dict[name]

class ColinerRec(object):
    def __init__(self, group_id, pair_id, name1, name2, score):
        self.group_id = group_id
        self.pair_id = pair_id
        self.name1 = name1
        self.name2 = name2
        self.score = score
        self.uniq_cds1 = None
        self.uniq_cds2 = None
        self.uniq_locus1 = None
        self.uniq_locus2 = None

    def find_uniq_id(self, genome1, genome2):
        self.uniq_cds1 = genome1.uniq_cds_map[self.name1]
        self.uniq_cds2 = genome2.uniq_cds_map[self.name2]
        self.uniq_locus1 = genome1.uniq_locus_map[self.name1]
        self.uniq_locus2 = genome2.uniq_locus_map[self.name2]

    @property
    def uniq_cds_pair(self):
        return (self.uniq_cds1, self.uniq_cds2)

    @property
    def uniq_locus_pair(self):
        return (self.uniq_locus1, self.uniq_locus2)


class Coliner(object):
    def __init__(self, genome1, genome2):
        self.genome1 = genome1
        self.genome2 = genome2
        self.coliner_dict = defaultdict(list)
        self.pair_dict = dict()

    def load_coliner(self, f_coliner):
        group_pattern = re.compile("\s*(\d+)-\s*(\d+)")
        for line in iterline(f_coliner):
            data = line.rstrip().split("\t")
            g_data, cds1, cds2, score = data
            group_id, pair_id = re.findall(group_pattern, g_data)[0]
            score = score.strip()
            cds1 = cds1.rstrip("_CDS")
            cds2 = cds2.rstrip("_CDS")
            if self.genome1.is_in_genome(cds1):
                if self.genome2.is_in_genome(cds2):
                    name1 = cds1
                    name2 = cds2
                else:
                    continue
            else:
                assert self.genome2.is_in_genome(cds1)
                if self.genome1.is_in_genome(cds2):
                    name1 = cds2
                    name2 = cds1
                else:
                    continue
            coliner_rec = ColinerRec(int(group_id), int(pair_id), name1, name2, score)
            coliner_rec.find_uniq_id(self.genome1, self.genome2)
            self.coliner_dict[coliner_rec.group_id].append(coliner_rec)
            self.pair_dict[(coliner_rec.name1, coliner_rec.name2)] = coliner_rec


class BlastRec(object):
    def __init__(self, name1, name2):
        self.name1 = name1
        self.name2 = name2
        self.uniq_cds1 = None
        self.uniq_cds2 = None
        self.uniq_locus1 = None
        self.uniq_locus2 = None
        self.fwd_match = False
        self.rev_match = False

    def find_uniq_id(self, genome1, genome2):
        self.uniq_cds1 = genome1.uniq_cds_map[self.name1]
        self.uniq_cds2 = genome2.uniq_cds_map[self.name2]
        self.uniq_locus1 = genome1.uniq_locus_map[self.name1]
        self.uniq_locus2 = genome2.uniq_locus_map[self.name2]

    def set_pos(self, fwd):
        if fwd:
            self.fwd_match = True
        else:
            self.rev_match = True

    @property
    def is_pair(self):
        if self.fwd_match and self.rev_match:
            return True
        return False

    @property
    def uniq_cds_pair(self):
        return (self.uniq_cds1, self.uniq_cds2)

    @property
    def uniq_locus_pair(self):
        return (self.uniq_locus1, self.uniq_locus2)

class Blast(object):
    def __init__(self, genome1, genome2):
        self.genome1 = genome1
        self.genome2 = genome2
        self.pair_dict = dict()

    def load_blast(self, f_blast):
        for line in iterline(f_blast):
            data = line.rstrip().split("\t")
            cds1 = data[0].rstrip("_CDS")
            cds2 = data[1].rstrip("_CDS")
            if self.genome1.is_in_genome(cds1):
                if self.genome2.is_in_genome(cds2):
                    name1 = cds1
                    name2 = cds2
                    fwd = True
                else:
                    continue
            else:
                assert self.genome2.is_in_genome(cds1)
                if self.genome1.is_in_genome(cds2):
                    name1 = cds2
                    name2 = cds1
                    fwd = False
                else:
                    continue
            key = (name1, name2)
            if key not in self.pair_dict:
                blast_rec = BlastRec(name1, name2)
                blast_rec.find_uniq_id(self.genome1, self.genome2)
                self.pair_dict[key] = blast_rec
            self.pair_dict[key].set_pos(fwd)


class HomologRec(object):
    TAG = {1: "SingleBlast", 2: "PairedBlast", 3: "Coliner"}

    def __init__(self, name):
        self.name = name
        self.level = 0

    def set_coliner(self):
        self.level = 3

    def set_pair_blast(self):
        self.level = max(self.level, 2)

    def set_single_blast(self):
        self.level = max(self.level, 1)

    def to_tuple(self):
        assert self.level > 0
        return (*self.name, self.TAG[self.level])


def scan_homolog_uniq_info(coliner_obj, blast_obj):
    cds_homolog = dict()
    locus_homolog = dict()

    for rec in coliner_obj.pair_dict.values():
        tmp_cds = rec.uniq_cds_pair
        if tmp_cds not in cds_homolog.keys():
            cds_homolog[tmp_cds] = HomologRec(tmp_cds)
        tmp_locus = rec.uniq_locus_pair
        if tmp_locus not in locus_homolog.keys():
            locus_homolog[tmp_locus] = HomologRec(tmp_locus)
        cds_homolog[tmp_cds].set_coliner()
        locus_homolog[tmp_locus].set_coliner()

    for rec in blast_obj.pair_dict.values():
        tmp_cds = rec.uniq_cds_pair
        if tmp_cds not in cds_homolog.keys():
            cds_homolog[tmp_cds] = HomologRec(tmp_cds)
        tmp_locus = rec.uniq_locus_pair
        if tmp_locus not in locus_homolog.keys():
            locus_homolog[tmp_locus] = HomologRec(tmp_locus)
        if rec.is_pair:
            cds_homolog[tmp_cds].set_pair_blast()
            locus_homolog[tmp_locus].set_pair_blast()
        else:
            cds_homolog[tmp_cds].set_single_blast()
            locus_homolog[tmp_locus].set_single_blast()

    cds_homolog_li = [rec.to_tuple() for _, rec in sorted(cds_homolog.items())]
    locus_homolog_li = [rec.to_tuple() for _, rec in sorted(locus_homolog.items())]
    return cds_homolog_li, locus_homolog_li


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--ref1", type=str, dest="ref1", metavar="ref1.cds.bed12", required=True)
    base_group.add_argument("--ref2", type=str, dest="ref2", metavar="ref2.cds.bed12", required=True)
    base_group.add_argument("--uniq-cds1", type=str, dest="uniq_cds1", metavar="uniq_cds1.tsv", required=True)
    base_group.add_argument("--uniq-cds2", type=str, dest="uniq_cds2", metavar="uniq_cds2.tsv", required=True)
    base_group.add_argument("--uniq-locus1", type=str, dest="uniq_locus1", metavar="uniq_locus1.tsv", required=True)
    base_group.add_argument("--uniq-locus2", type=str, dest="uniq_locus2", metavar="uniq_locus2.tsv", required=True)
    base_group.add_argument("--assemble1", type=str, dest="assemble1", metavar="assemble1", required=True)
    base_group.add_argument("--assemble2", type=str, dest="assemble2", metavar="assemble2", required=True)
    base_group.add_argument("--coliner", type=str, dest="coliner", metavar="coliner.txt", required=True)
    base_group.add_argument("--blast", type=str, dest="blast", metavar="blast.tsv", required=True)
    base_group.add_argument("--out-homolog-cds", type=str, dest="out_homolog_cds", metavar="homolog.cds.tsv",
                            required=True)
    base_group.add_argument("--out-homolog-locus", type=str, dest="out_homolog_locus", metavar="homolog.locus.tsv",
                            required=True)
    base_group.add_argument("--out-coliner-pos", type=str, dest="out_coliner_pos", metavar="coliner.pos.tsv",
                            required=True)
    base_group.add_argument("--out-coliner-locus", type=str, dest="out_coliner_locus", metavar="coliner.locus.pair.tsv",
                            required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_ref1 = args.ref1
    f_ref2 = args.ref2
    f_uniq_cds1 = args.uniq_cds1
    f_uniq_cds2 = args.uniq_cds2
    f_uniq_locus1 = args.uniq_locus1
    f_uniq_locus2 = args.uniq_locus2
    assemble1 = args.assemble1
    assemble2 = args.assemble2
    f_coliner = args.coliner
    f_blast = args.blast
    f_out_homolog_cds = args.out_homolog_cds
    f_out_homolog_locus = args.out_homolog_locus
    f_out_coliner_pos = args.out_coliner_pos
    f_out_coliner_locus = args.out_coliner_locus

    genome1 = GenomeInfo(assemble1, f_ref1, f_uniq_cds1, f_uniq_locus1)
    genome2 = GenomeInfo(assemble2, f_ref2, f_uniq_cds2, f_uniq_locus2)

    coliner_obj = Coliner(genome1, genome2)
    coliner_obj.load_coliner(f_coliner)

    f_locus = open(f_out_coliner_locus, "w")
    out_coliner_locus_header = "ColinerID\tPairID\tAssemble1\tUniqLocus1\tAssemble2\tUniqLocus2\tScore\n"
    f_locus.write(out_coliner_locus_header)
    f_pos = open(f_out_coliner_pos, "w")
    out_coliner_pos_header = "ColinerID\tAssemble1\tChrom1\tStart1\tEnd1\tAssemble2\tChrom2\tStart2\tEnd2\n"
    f_pos.write(out_coliner_pos_header)

    for coliner_indx, (group_id, rec_li) in enumerate(sorted(coliner_obj.coliner_dict.items())):
        first_rec = rec_li[0]
        last_rec = rec_li[-1]
        first_cds1 = genome1.get_cds(first_rec.name1)
        first_cds2 = genome2.get_cds(first_rec.name2)
        last_cds1 = genome1.get_cds(last_rec.name1)
        last_cds2 = genome2.get_cds(last_rec.name2)

        assert first_cds1.chrom == last_cds1.chrom, (first_cds1.name, last_cds1.name)
        assert first_cds2.chrom == last_cds2.chrom, (first_cds2.name, last_cds2.name)

        data = [
            coliner_indx + 1,
            assemble1, first_cds1.chrom, min(first_cds1.chromStart, last_cds1.chromStart),
            max(first_cds1.chromEnd, last_cds1.chromEnd), 
            assemble2, first_cds2.chrom, min(first_cds2.chromStart, last_cds2.chromStart),
            max(first_cds2.chromEnd, last_cds2.chromEnd)
        ]
        f_pos.write("\t".join(list(map(str, data))) + "\n")

        for pair_indx, rec in enumerate(rec_li):
            data = [coliner_indx + 1, pair_indx + 1, assemble1, rec.uniq_locus1, assemble2, rec.uniq_locus2, rec.score]
            f_locus.write("\t".join(list(map(str, data))) + "\n")

    blast_obj = Blast(genome1, genome2)
    blast_obj.load_blast(f_blast)

    cds_homolog_li, locus_homolog_li = scan_homolog_uniq_info(coliner_obj, blast_obj)

    with open(f_out_homolog_cds, "w") as f:
        out_homolog_cds_header = "Assemble1\tUniqCDS1\tAssemble2\tUniqCDS2\tType\n"
        f.write(out_homolog_cds_header)
        for name1, name2, tag in cds_homolog_li:
            data = [assemble1, name1, assemble2, name2, tag]
            f.write("\t".join(data) + "\n")

    with open(f_out_homolog_locus, "w") as f:
        out_homolog_locus_header = "Species1\tUniqLocus1\tSpecies2\tUniqLocus2\tType\n"
        f.write(out_homolog_locus_header)
        for name1, name2, tag in locus_homolog_li:
            data = [assemble1, name1, assemble2, name2, tag]
            f.write("\t".join(data) + "\n")
            

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
