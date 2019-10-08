import os
curDir = os.getcwd()
MSSEARCH_BASE = os.path.join(curDir, "results", "MSSearch")
TOOLS_BASE = os.path.join(curDir, "tools")


DATABASE = {
    "gosHir":[
        "gosHir.BGI_1_0", "gosHir.JGI_1_1"
    ],
    "theCac":[
        "theCac.cocoa_v2"
    ]
}
SAMPLES = {
    "gosHir":[
        "P180076_F2", "P180076_F4", "P180076_F5"
    ], 
    "theCac":[
        "Fraction_after_buffer_extraction", "Salt_soluble_fraction",
        "Unfractionated_sample", "Water_soluble_fraction"
    ]
}




MAKE_DECOY_DATABASE_BASE = os.path.join(MSSEARCH_BASE, "make_decoy_database")
rule make_decoy_database:
    input:
        faa = os.path.join(curDir, "data", "annotation", "{database}.cds.faa"),
    output:
        custom_faa = os.path.join(MAKE_DECOY_DATABASE_BASE, "{database}.cds.custom.faa"),
        customDe_fa = os.path.join(MAKE_DECOY_DATABASE_BASE, "{database}.cds.custom_concatenated_target_decoy.fasta"),
    params:
        make_custom_faa_py = os.path.join(curDir, "scripts", "MSSearch", "makeCustomFaa.py"),
        searchGUI = os.path.join(TOOLS_BASE, "searchGUI", "SearchGUI-3.3.15", "SearchGUI-3.3.15.jar"),
    shell:
        """
        python {params.make_custom_faa_py} -i {input.faa} -o {output.custom_faa}
        java -Djava.awt.headless=true -cp {params.searchGUI} eu.isas.searchgui.cmd.FastaCLI -in {output.custom_faa} -decoy
        """


MAKE_PARAMS_BASE = os.path.join(MSSEARCH_BASE, "make_params")
rule make_params:
    input:
        customDe_fa = rules.make_decoy_database.output.customDe_fa,
    output:
        params = os.path.join(MAKE_PARAMS_BASE, "{database}.params.par"),
    params:
        searchGUI = os.path.join(TOOLS_BASE, "searchGUI", "SearchGUI-3.3.15", "SearchGUI-3.3.15.jar"),
        out_path = os.path.join(MAKE_PARAMS_BASE, "{database}.params"),
    shell:
        """
        java -Djava.awt.headless=true -cp {params.searchGUI} eu.isas.searchgui.cmd.IdentificationParametersCLI -out {params.out_path} -db {input.customDe_fa}
        """

RUN_SEARCHGUI_BASE = os.path.join(MSSEARCH_BASE, "run_searchGUI")
rule run_searchGUI:
    input:
        params = rules.make_params.output.params,
        spectrum_files = os.path.join(curDir, "results", "ConvertMGF", "mgf", "{sample}.mgf"),
    output:
        xtandem_out = os.path.join(RUN_SEARCHGUI_BASE, "{species}", "{database}", "{sample}", "searchgui_out.zip"),
    params:
        searchGUI = os.path.join(TOOLS_BASE, "searchGUI", "SearchGUI-3.3.15", "SearchGUI-3.3.15.jar"),
        output_folder = os.path.join(RUN_SEARCHGUI_BASE, "{species}", "{database}", "{sample}"),
        temp_folder = os.path.join(RUN_SEARCHGUI_BASE, "{species}", "{database}", "{sample}", "run_searchGUI.tmp"),
    log:
        os.path.join(RUN_SEARCHGUI_BASE, "{species}", "{database}", "{sample}", "run_searchGUI.log"),
    threads: 4
    shell:
        """
        java -Djava.awt.headless=true -Xmx100G -cp {params.searchGUI} eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.spectrum_files} -output_folder {params.output_folder} -id_params {input.params} -xtandem 1 -threads {threads} -log {log} -temp_folder {params.temp_folder}
        """

RUN_PEPTIDESHAKER_BASE = os.path.join(MSSEARCH_BASE, "run_peptideShaker")
rule run_peptideShaker:
    input:
        params = rules.make_params.output.params,
        xtandem_out = rules.run_searchGUI.output.xtandem_out,
        spectrum_files = rules.run_searchGUI.input.spectrum_files,
    output:
        peptideShaker_out = os.path.join(RUN_PEPTIDESHAKER_BASE, "{species}", "{database}", "{sample}", "peptideShaker.cpsx"),
    params:
        peptideShaker = os.path.join(TOOLS_BASE, "peptideShaker", "PeptideShaker-1.16.40", "PeptideShaker-1.16.40.jar"),
        sample_name = "{sample}",
        temp_folder = os.path.join(RUN_PEPTIDESHAKER_BASE, "{species}", "{database}", "{sample}", "run_peptideShaker.tmp"),
    log:
        os.path.join(RUN_PEPTIDESHAKER_BASE, "{species}", "{database}", "{sample}", "run_peptideShaker.log"),
    threads: 4
    shell:
        """
        java -Djava.awt.headless=true -Xmx100G -cp {params.peptideShaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment myExperiment -sample {params.sample_name} -replicate 1 -identification_files {input.xtandem_out} -spectrum_files {input.spectrum_files} -id_params {input.params} -out {output.peptideShaker_out} -threads {threads} -log {log} -temp_folder {params.temp_folder}
        """

RUN_REPORTCLI_BASE = os.path.join(MSSEARCH_BASE, "run_reportCLI")
rule run_reportCLI:
    input:
        peptideShaker_out = rules.run_peptideShaker.output.peptideShaker_out
    output:
        reportCLI_out = os.path.join(RUN_REPORTCLI_BASE, "{species}", "{database}", "{sample}", "ReportCLI/"),
        pep_txt = os.path.join(RUN_REPORTCLI_BASE, "{species}", "{database}", "{sample}", "ReportCLI", "myExperiment_{sample}_1_Default_Peptide_Report.txt"),
    params:
        peptideShaker = os.path.join(TOOLS_BASE, "peptideShaker", "PeptideShaker-1.16.40", "PeptideShaker-1.16.40.jar"),
        temp_folder = os.path.join(RUN_REPORTCLI_BASE, "{species}", "{database}", "{sample}", "run_reportCLI.tmp"),
    log:
        os.path.join(RUN_REPORTCLI_BASE, "{species}", "{database}", "{sample}", "run_reportCLI.log"),
    threads: 28
    shell:
        """
        java -Djava.awt.headless=true -Xmx150G -cp {params.peptideShaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptideShaker_out} -out_reports {output.reportCLI_out} -reports 0,1,3,6,9 -temp_folder {params.temp_folder} 1>{log} 2>&1
        """

PEP_TO_BED_BASE = os.path.join(MSSEARCH_BASE, "pep_to_bed")
rule pep_to_bed:
    input:
        rules.run_reportCLI.output.reportCLI_out,
        pep_txt = rules.run_reportCLI.output.pep_txt,
        db_bed12 = os.path.join(curDir, "data", "annotation", "{database}.cds.bed12"),
        db_faa = os.path.join(curDir, "data", "annotation", "{database}.cds.faa"),
        db_genome = os.path.join(curDir, "data", "annotation", "{database}.genome.fa"),
    output:
        out_bed12 = os.path.join(PEP_TO_BED_BASE, "{species}", "{database}", "{sample}.bed12"),
    params:
        py = os.path.join(curDir, "scripts", "MSSearch", "pepToBed.py"),
    shell:
        """
        python {params.py} -i {input.pep_txt} -b {input.db_bed12} -f {input.db_faa} -g {input.db_genome} -o {output.out_bed12}
        """


BED_TO_TSV_BASE = os.path.join(MSSEARCH_BASE, "bed_to_tsv")
rule bed_to_tsv:
    input:
        in_bed12 = rules.pep_to_bed.output.out_bed12,
    output:
        out_tsv = os.path.join(BED_TO_TSV_BASE, "{species}", "{database}", "{sample}.tsv"),
    params:
        py = os.path.join(curDir, "scripts", "MSSearch", "bedToTsv.py"),
    shell:
        """
        python {params.py} -i {input.in_bed12} -o {output.out_tsv}
        """

def expand_database(ruleStr):
    res = list()
    for tmp_species in DATABASE:
        for tmp_database in DATABASE[tmp_species]:
            res.append(ruleStr.format(database = tmp_database))
    return(res)


def expand_sample(ruleStr):
    res = list()
    for tmp_species in DATABASE:
        for tmp_database in DATABASE[tmp_species]:
            for tmp_sample in SAMPLES[tmp_species]:
                res.append(ruleStr.format(species = tmp_species, database = tmp_database, sample = tmp_sample))
    return(res)

rule MSSearch:
    input:
        customDe_fa = expand_database(rules.make_decoy_database.output.customDe_fa),
        params = expand_database(rules.make_params.output.params),
        xtandem_out = expand_sample(rules.run_searchGUI.output.xtandem_out),
        peptideShaker_out = expand_sample(rules.run_peptideShaker.output.peptideShaker_out),
        reportCLI_out = expand_sample(rules.run_reportCLI.output.reportCLI_out),
        out_bed12 = expand_sample(rules.pep_to_bed.output.out_bed12),
        out_tsv = expand_sample(rules.bed_to_tsv.output.out_tsv),
