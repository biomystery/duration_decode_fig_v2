while read g; do
  echo $g
  Rscript fit_rule_shutdown_all.R $g > "./3rd_fit_9pars_measured_deg/"$g".output.txt" &
done < ./data/glist_final.csv
