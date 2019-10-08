source activate py35
python3 ../../summary_blast_out.py -i input.blast.out -o res.tsv
if [[ -n $(diff exp_res.tsv res.tsv) ]]; then
    diff exp_res.tsv res.tsv
    exit 1
fi
rm res.tsv
