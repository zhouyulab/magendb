#!/usr/bin/env python3

import os, sys, argparse
from Bio import SeqIO

def recordBed12(f_bed12):
    bed12_ref = {}
    with open(f_bed12, "r") as f:
        for line_data in f:
            line = line_data.strip("\n").split("\t")
            trans_id = line[3]
            bed12_ref[trans_id] = line
    return(bed12_ref)

def recordFasta(f_fasta):
    pro_ref = {}
    with open(f_fasta, "r") as f:
        for r in SeqIO.parse(f, "fasta"):
            pro_ref[r.id] = str(r.seq)
    return(pro_ref)

def getIndex(site, blockSizeList, trans_strand):
    #site = 4758
    #trans_strand = "-"
    if trans_strand == "+":
        exon_sum_init = 0
        for i in range(len(blockSizeList)):
            exon_sum = exon_sum_init + blockSizeList[i]
            if site >= exon_sum_init and site <= exon_sum:
                index = i
                break
            exon_sum_init = exon_sum
        exon_residue = site - sum(blockSizeList[:i])
        return(index, exon_residue)
    elif trans_strand == "-":
        blockSizeList = blockSizeList[::-1]
        exon_sum_init = 0
        for i in range(len(blockSizeList)):
            exon_sum = exon_sum_init + blockSizeList[i]
            if site >= exon_sum_init and site <= exon_sum:
                index = len(blockSizeList) - i - 1
                break
            exon_sum_init = exon_sum
        exon_residue = site - sum(blockSizeList[:i])
        return(index, exon_residue)
    


def makePepBed12(trans_info, start_rel_index, start_exon_residue, end_rel_index, end_exon_residue, trans_strand):
    ### make transcript bed12
    transcript_bed12 = "%s/transcript.bed12"%tmp_dir
    handle = open(transcript_bed12, "w")
    handle.write("\t".join(trans_info) + "\n")
    handle.close()

    ### transcript bed12 to bed6
    transcript_bed6 = "%s/transcript.bed6"%tmp_dir
    command = "bedtools bed12tobed6 -i %s > %s"%(transcript_bed12, transcript_bed6)
    os.system(command)

    ### make peptide bed6
    pep_bed6 = "%s/pep.bed6"%tmp_dir
    handle_pep = open(pep_bed6, "w")
    with open(transcript_bed6, "r") as f:
        transcript_exons = f.readlines()

    if start_rel_index == end_rel_index:
        target_exon = transcript_exons[start_rel_index].strip("\n").split("\t")
        start_site, end_site = int(target_exon[1]), int(target_exon[2])
        if trans_strand == "+":
            target_exon[1] = str(start_site + start_exon_residue)
            target_exon[2] = str(start_site + end_exon_residue)
        elif trans_strand == "-":
            target_exon[1] = str(end_site - end_exon_residue)
            target_exon[2] = str(end_site - start_exon_residue)
        handle_pep.write("\t".join(target_exon) + "\n")
    elif start_rel_index != end_rel_index:
        start_exon = transcript_exons[start_rel_index].strip("\n").split("\t")
        end_exon = transcript_exons[end_rel_index].strip("\n").split("\t")
        if trans_strand == "+":
            start_exon[1] = str(int(start_exon[1]) + start_exon_residue)
            end_exon[2] = str(int(end_exon[1]) + end_exon_residue)
            other_exon = transcript_exons[start_rel_index+1:end_rel_index]
        elif trans_strand == "-":
            start_exon[2] = str(int(start_exon[2]) - start_exon_residue)
            end_exon[1] = str(int(end_exon[2]) - end_exon_residue)
            other_exon = transcript_exons[end_rel_index+1:start_rel_index]
        handle_pep.write("\t".join(start_exon) + "\n")
        handle_pep.write("\t".join(end_exon) + "\n")
        for o in range(len(other_exon)):
            handle_pep.write(other_exon[o])
    handle_pep.close()

    ### bed6 to bed12
    pep_gtf = "%s/pep.gtf"%tmp_dir
    handle_gtf = open(pep_gtf, "w")
    with open(pep_bed6, "r") as f:
        for line_data in f:
            line = line_data.strip("\n").split("\t")
            if line[1] == line[2]:
                continue
            exon_info = [line[0], "peptide", "exon", str(int(line[1])+1), line[2], "0", line[5], ".", 'gene_id "peptide"; transcript_id "%s";'%line[3]]
            handle_gtf.write("\t".join(exon_info) + "\n")
    handle_gtf.close()

    pep_genePred = "%s/pep.genePred"%tmp_dir
    command_gtfToGenepred = "gtfToGenePred %s %s"%(pep_gtf, pep_genePred)
    os.system(command_gtfToGenepred)
        
    pep_bed12 = "%s/pep.bed12"%tmp_dir
    command_genePredToBed12 = "genePredToBed %s %s"%(pep_genePred, pep_bed12)
    os.system(command_genePredToBed12)
    return(pep_bed12)
    

class peptide(object):
    def __init__(self, line_data): 
        self.line = line_data.strip("\n").split("\t")
        self.index = self.line[0]
        self.proteins = self.line[1].split(";")
        self.sequence = self.line[5]
        self.modification = self.line[6]
        self.PSMs = self.line[13]
        self.validation = self.line[16]

    def getRelPosition_pro(self, pro_ref, pro):
        pro_seq = pro_ref[pro]
        rel_p, i = [], 0
        while pro_seq.find(self.sequence, i) != -1:
            
            rel_p.append(pro_seq.find(self.sequence, i)+1)
            i = pro_seq.find(self.sequence, i)+1 
        return(rel_p)        

    def getRelPosition_nuc(self, p): 
        start_rel = (p - 1)*3
        end_rel = (p -1 + len(self.sequence))*3
        return(start_rel, end_rel)
       

    def makeBed12(self, trans_info, start_rel, end_rel):
        blockSizeList = list(map(int, trans_info[10].split(",")))
        trans_strand = trans_info[5]

        ### get start_rel_index and end_rel_index
        start_rel_index, start_exon_residue = getIndex(start_rel, blockSizeList, trans_strand)
        end_rel_index, end_exon_residue  = getIndex(end_rel, blockSizeList, trans_strand)
        #if start_rel_index+1 < end_rel_index:
        #    print(start_rel_index, start_exon_residue, end_rel_index, end_exon_residue, trans_info)
        #    print(self.index)
   
        ### make peptide bed12
        pep_bed12 = makePepBed12(trans_info, start_rel_index, start_exon_residue, end_rel_index, end_exon_residue, trans_strand)
        return(pep_bed12)

    def recoverPep(self, pep_bed12, f_genome):
        pep_fa = "%s/pep.fa"%tmp_dir
        command_getfasta = "bedtools getfasta -fi %s -s -split -bed %s -fo %s"%(f_genome, pep_bed12, pep_fa)
        os.system(command_getfasta)
        
        with open(pep_fa, "r") as f:
            r = SeqIO.read(f, "fasta")
            pep_seq = r.seq.translate()
        return(pep_seq) 

    def output(self, pep_bed12, pep_seq, p, handle_out):
         with open(pep_bed12, "r") as f:
             for line_data in f:
                 line = line_data.strip("\n").split("\t")
                 modify = self.modification.split("-")
                 pep_seq = "%s-%s-%s"%(modify[0], pep_seq, modify[2])
                 line[3] = "%s|%s|position:%s"%(line[3], pep_seq, p)
                 line[4] = self.PSMs
                 handle_out.write("\t".join(line) + "\n")


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.txt", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.bed12", required=True)
    base_group.add_argument("-b", "--bed12", type=str, dest="bed12", metavar="database.bed12", required=True)
    base_group.add_argument("-f", "--fasta", type=str, dest="fasta", metavar="database.faa", required=True)
    base_group.add_argument("-g", "--genome", type=str, dest="genome", metavar="genome.fa", required=True)

    return(parser.parse_args(args))


def main():
    args = sys.argv[1:]

    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_bed12 = args.bed12
    f_fasta = args.fasta
    f_genome = args.genome

    ### create tmp_dir
    global tmp_dir
    tmp_dir = "%s.tmp"%f_out
    if os.path.exists(tmp_dir) == True:
        os.system("rm -rf %s"%tmp_dir)
    os.system("mkdir %s"%tmp_dir)

    ### record all transcripts bed12
    bed12_ref = recordBed12(f_bed12)

    ### record all protein fasta
    pro_ref = recordFasta(f_fasta)

    ### 
    handle_out = open(f_out, "w")

    ###
    with open(f_in, "r") as f:
        next(f)
        for line_data in f:
            #print(line_data)
            ### definite a peptide
            pep = peptide(line_data)
 
            if pep.validation != "Confident":
                continue

            for pro in pep.proteins:
                pro = pro.replace(" ", "")
                ### get peptide relative position in protein
                rel_p = pep.getRelPosition_pro(pro_ref, pro)

                ### get transcript bed12 information
                pep_transID = pro
                trans_info = bed12_ref[pep_transID]

                for p in rel_p:
                    ### get peptide relative position in transcript
                    start_rel, end_rel = pep.getRelPosition_nuc(p)

                    ### make peptide bed12
                    pep_bed12 = pep.makeBed12(trans_info, start_rel, end_rel)

                    ### recover peptide sequence
                    pep_seq = pep.recoverPep(pep_bed12, f_genome)
                
                    ### output
                    pep.output(pep_bed12, pep_seq, p, handle_out)
            #break
    os.system("rm -rf %s"%tmp_dir)
    handle_out.close()
   


if __name__ == "__main__":
    main()
