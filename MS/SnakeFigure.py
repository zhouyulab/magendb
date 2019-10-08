import os
curDir = os.getcwd()
FIGURE_BASE = os.path.join(curDir, "results", "figure")

include: "SnakeMSSearch.py"

rule data_summary:
    input:
        tsv = rules.MSSearch.input.out_tsv,
    output:
        tsv_gene = os.path.join(FIGURE_BASE, "data_summary_gene.tsv"),
        tsv_pep = os.path.join(FIGURE_BASE, "data_summary_peptide.tsv"),
    shell:
        """
        python scripts/figure/data_summary.py {input.tsv} {output.tsv_gene} {output.tsv_pep}
        """

rule pep_per_gene:
    input:
        tsv = rules.data_summary.output.tsv_gene, 
    output:
        pdf = os.path.join(FIGURE_BASE, "pep_per_gene.pdf"),
        tsv = os.path.join(FIGURE_BASE, "pep_per_gene.tsv"),
    shell:
        """
        Rscript scripts/figure/pep_per_gene.R {input.tsv} {output.pdf} {output.tsv}
        """  

rule PSM_per_gene:
    input:
        tsv = rules.data_summary.output.tsv_gene,
    output:
        pdf = os.path.join(FIGURE_BASE, "PSM_per_gene.pdf"),
        tsv = os.path.join(FIGURE_BASE, "PSM_per_gene.tsv"),
    shell:
        """
        Rscript scripts/figure/PSM_per_gene.R {input.tsv} {output.pdf} {output.tsv}
        """


rule Figure:
    input:
        rules.data_summary.output.tsv_gene,
        rules.pep_per_gene.output.pdf,
        rules.PSM_per_gene.output.pdf,
