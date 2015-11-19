## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout

mkdir $ROUT

### Run R scripts


for i in 'drosophila_node_nonull'
do 
  for j in 'kallisto' 'htseq'
  do

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 count_method='${j}' simulation='${i}' dispersion_common=TRUE results_common=TRUE disp_mode_list=c('grid','grid') disp_moderation_list=c('none','common')" $RCODE/sim5_drimseq_0_3_1_run.R $ROUT/sim5_drimseq_0_3_1_run_${i}_${j}_grid.Rout
    
    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 count_method='${j}' simulation='${i}' dispersion_common=FALSE results_common=FALSE disp_mode_list=c('optimize','optim','constrOptim') disp_moderation_list=c('none','none','none')" $RCODE/sim5_drimseq_0_3_1_run.R $ROUT/sim5_drimseq_0_3_1_run_${i}_${j}_optim.Rout

  done
done


for i in 'drosophila_node_nonull'
do 
  for j in 'htseq_prefiltered15'
  do

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 count_method='${j}' simulation='${i}' dispersion_common=TRUE results_common=TRUE disp_mode_list=c('grid','grid') disp_moderation_list=c('none','common')" $RCODE/sim5_drimseq_0_3_1_run.R $ROUT/sim5_drimseq_0_3_1_run_${i}_${j}_grid.Rout
    
  done
done



for i in 'hsapiens_node_nonull'
do 
  for j in 'kallisto' 'htseq'
  do

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 count_method='${j}' simulation='${i}' dispersion_common=TRUE results_common=TRUE disp_mode_list=c('grid','grid') disp_moderation_list=c('none','common')" $RCODE/sim5_drimseq_0_3_1_run.R $ROUT/sim5_drimseq_0_3_1_run_${i}_${j}_grid.Rout
    
    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 count_method='${j}' simulation='${i}' dispersion_common=FALSE results_common=FALSE disp_mode_list=c('optimize','optim','constrOptim') disp_moderation_list=c('none','none','none')" $RCODE/sim5_drimseq_0_3_1_run.R $ROUT/sim5_drimseq_0_3_1_run_${i}_${j}_optim.Rout

  done
done


for i in 'hsapiens_node_nonull'
do 
  for j in 'htseq_prefiltered15'
  do

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 count_method='${j}' simulation='${i}' dispersion_common=TRUE results_common=TRUE disp_mode_list=c('grid','grid') disp_moderation_list=c('none','common')" $RCODE/sim5_drimseq_0_3_1_run.R $ROUT/sim5_drimseq_0_3_1_run_${i}_${j}_grid.Rout

  done
done



for i in 'hsapiens_withde_nonull'
do 
  for j in 'kallisto' 'htseq'
  do

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 count_method='${j}' simulation='${i}' dispersion_common=TRUE results_common=TRUE disp_mode_list=c('grid','grid') disp_moderation_list=c('none','common')" $RCODE/sim5_drimseq_0_3_1_run.R $ROUT/sim5_drimseq_0_3_1_run_${i}_${j}_grid.Rout

  done
done







