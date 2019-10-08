source activate py35
python3 ../../summary_utr_scan_rg4_qgrs.py -i rg4.out -r iso.bed12 -o res.bed12
if [[ -n $(diff exp_res.bed12 res.bed12) ]]; then
    diff exp_res.bed12 res.bed12
    exit 1
fi
rm res.bed12
