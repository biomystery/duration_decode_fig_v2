awk -v FS=',' '(NR>1){print $1}' ../data/model_v2_pars.csv| while read g; do
  echo $g
  Rscript fit_rule_shutdown_all.R $g hill_1/ > "./hill_1/"$g".output.txt" &
done 
