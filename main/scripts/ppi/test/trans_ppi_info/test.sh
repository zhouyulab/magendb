source activate py35
python3 ../../trans_ppi_info.py --coliner colinear.tsv --ppi ppi.tsv --string-ref string_ref.tsv --plant-ref plant_ref.tsv -o res.tsv
if [[ -n $(diff exp_res.tsv res.tsv) ]]; then
    diff exp_res.tsv res.tsv
    exit 1
fi
rm res.tsv
