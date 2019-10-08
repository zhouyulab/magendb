source activate py35
python3 ../../parse_planttfdb_blast_res.py -i input.blast.out -o res.tsv
if [[ -n $(diff exp_res.tsv res.tsv) ]]; then
    diff exp_res.tsv res.tsv
    exit 1
fi
rm res.tsv
