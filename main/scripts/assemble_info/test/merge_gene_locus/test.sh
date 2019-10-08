source activate py35
python3 ../../merge_gene_locus.py -i A1.bed12 A2.bed12 -t A1 A2 -s test --gene-locus locus.bed12 --gene-locus-map locus_map.tsv
if [[ -n $(diff exp_locus.bed12 locus.bed12) ]]; then
    diff exp_locus.bed12 locus.bed12
    exit 1
fi
if [[ -n $(diff exp_locus_map.tsv locus_map.tsv) ]]; then
    diff exp_locus_map.tsv locus_map.tsv
    exit 1
fi
rm locus.bed12 locus_map.tsv
