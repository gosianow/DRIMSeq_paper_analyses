## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout

mkdir $ROUT

### Run R scripts


for i in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for j in 'kallisto' 'htseq' 'htseq_prefiltered15'
  do
    
    echo "${i}_${j}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${j}' simulation='${i}'" $RCODE/sim5_drimseq_0_3_1_comparison_run.R $ROUT/sim5_drimseq_0_3_1_comparison_run_${i}_${j}.Rout

  done
done



R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/sim5_drimseq_0_3_1_comparison_combined.R $ROUT/sim5_drimseq_0_3_1_comparison_combined.Rout


