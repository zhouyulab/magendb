#!/usr/bin/env python3

import os, sys, argparse

def cds_map_gene():
    map_info = dict()
    anno_base = os.path.join("data", "annotation")
    gene_locus_list = ["theCac.cocoa_v2.gene_locus_map.tsv", "gosHir.BGI_1_0.gene_locus_map.tsv", "gosHir.JGI_1_1.gene_locus_map.tsv"]
    for g in gene_locus_list:
        g_path = os.path.join(anno_base, g)
        with open(g_path, "r") as f:
            for line_data in f:
                line = line_data.strip("\n").split("\t")
                if line[0] == "GeneLocus":
                    continue
                
                GeneLocus, Tid = line[0], line[2]
                map_info[Tid] = GeneLocus
    return(map_info)


def main():
    file_list = sys.argv[1:11]
    f_out_1 = sys.argv[11]
    f_out_2 = sys.argv[12]

    ### map cds to gene
    map_info = cds_map_gene()
   
    ### count gene peptides and PSMs 
    gene = dict()
    peptide = dict()

    for file in file_list:
        with open(file, "r") as f:
            for line_data in f:
                line = line_data.strip("\n").split("\t")
                if line[0] == "gene_id":
                    continue

                ### gene
                gene_id = line[0].split("_CDS")[0]
                gene_locus = map_info[gene_id]
                pep_nums = line[1].split(";")
                PSMs = sum(list(map(int, line[3].split(";"))))

                if gene_locus not in gene:
                    gene[gene_locus] = [pep_nums, PSMs]
                else:
                    gene[gene_locus][0] += pep_nums
                    gene[gene_locus][1] += PSMs
     
                ### peptide
                pep_list = line[1].split(";") 
                PSMs = list(map(int, line[3].split(";")))
             
                if gene_locus not in peptide:
                    peptide[gene_locus] = dict()
                for i, p in enumerate(pep_list):
                    if p not in peptide[gene_locus]:
                        peptide[gene_locus][p] = PSMs[i]
                    else:
                        peptide[gene_locus][p] += PSMs[i]

    ### output
    with open(f_out_1, "w") as f:
        ### make header
        header = ["gene_locus", "peptide_per_gene", "PSMs_per_gene"]
        f.write("\t".join(header) + "\n")
        for gene_locus in gene:
            info = [gene_locus, len(set(gene[gene_locus][0])), gene[gene_locus][1]]
            info = list(map(str, info))

            f.write("\t".join(info) + "\n")

    with open(f_out_2, "w") as f:
        ### make header
        header = ["gene_locus", "peptide", "PSMs_total"]
        f.write("\t".join(header) + "\n")
        for gene_locus in peptide:
            for p in peptide[gene_locus]:
                info = [gene_locus, p, peptide[gene_locus][p]]
                info = list(map(str, info))

                f.write("\t".join(info) + "\n")    



if __name__ == "__main__":
    main()

