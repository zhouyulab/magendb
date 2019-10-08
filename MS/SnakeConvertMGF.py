import os
curDir = os.getcwd()

SPECIES = ["gosHir", "theCac"]
SAMPLES = {
    "gosHir":[
        "P180076_F1", "P180076_F2", "P180076_F3", "P180076_F4", "P180076_F5"
    ], 
    "theCac":[
        "Fraction_after_buffer_extraction", "Salt_soluble_fraction",
        "Unfractionated_sample", "Water_soluble_fraction"
    ]
}

CONVERTMGF_BASE = os.path.join(curDir, "results", "convertMGF")
rule msconvert:
    input:
        raw_file = os.path.join(curDir, "data", "raw", "{sample}.raw"),
    output:
        mgf_file = os.path.join(curDir, "results", "ConvertMGF", "mgf", "{sample}.mgf"),
    params:
        docker_path = curDir,
    threads:1
    shell:
        """
        docker run --rm -e WINEDEBUG=-all -v {params.docker_path}:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert --mgf /data/data/raw/{wildcards.sample}.raw --outfile /data/results/ConvertMGF/mgf/{wildcards.sample}.mgf -o /data/results/ConvertMGF/mgf/ ; chown wukai:zhoulab {output.mgf_file}
        """

def expand_convertMGF(ruleStr):
    res = list()
    for tmp_species in SPECIES:
        for tmp_sample in SAMPLES[tmp_species]:
            res.append(ruleStr.format(sample = tmp_sample))
    return(res)

rule ConvertMGF:
    input:
        expand_convertMGF(rules.msconvert.output.mgf_file),
