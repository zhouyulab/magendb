source activate py35
python3 ../../summary_TMHMM_results.py -i input.out -r cds.bed12 -o res.bed12
if [[ -n $(diff exp_res.bed12 res.bed12) ]]; then
    diff exp_res.bed12 res.bed12
    exit 1
fi
rm res.bed12
