## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_brooks_pasilla
RWD=/home/Shared/data/seq/brooks_pasilla
ROUT=$RWD/Rout


# mkdir $ROUT

## Run R scripts


for model in 'model_full' 'model_full_glm' 'model_full_paired'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${model}_${count_method}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' model='${model}'" $RCODE/brooks_positive_controls.R $ROUT/brooks_positive_controls_run_${model}_${count_method}.Rout

  done
done


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/brooks_positive_controls_summary.R $ROUT/brooks_positive_controls_summary.Rout


