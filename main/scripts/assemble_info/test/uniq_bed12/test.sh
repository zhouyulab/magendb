source activate py35
python3 ../../uniq_bed12.py -i A1.bed12 A2.bed12 -t A1 A2 --uniq uniq.bed12 --uniq-map uniq_map.tsv
if [[ -n $(diff exp_uniq.bed12 uniq.bed12) ]]; then
    diff exp_uniq.bed12 uniq.bed12
    exit 1
fi
if [[ -n $(diff exp_uniq_map.tsv uniq_map.tsv) ]]; then
    diff exp_uniq_map.tsv uniq_map.tsv
    exit 1
fi
rm uniq.bed12 uniq_map.tsv
