source activate py35
python3 ../../merge_ppi_score.py -i A1.score.tsv A2.score.tsv --assemble A1 A2 --uniq-locus-map locus.map.tsv --out-prot prot.ppi.tsv --out-locus locus.ppi.tsv
if [[ -n $(diff exp_prot.ppi.tsv prot.ppi.tsv) ]]; then
    diff exp_prot.ppi.tsv prot.ppi.tsv
    exit 1
fi
if [[ -n $(diff exp_locus.ppi.tsv locus.ppi.tsv) ]]; then
    diff exp_locus.ppi.tsv locus.ppi.tsv
    exit 1
fi

rm prot.ppi.tsv locus.ppi.tsv
