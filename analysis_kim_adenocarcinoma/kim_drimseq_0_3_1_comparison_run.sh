## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout

# mkdir $ROUT

## Run R scripts


for model in 'model_full' 'model_null_normal1' 'model_null_normal2' 'model_null_tumor1' 'model_null_tumor2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${model}_${count_method}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' model='${model}'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run_${model}_${count_method}.Rout

  done
done



R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/kim_summary.R $ROUT/kim_summary.Rout


