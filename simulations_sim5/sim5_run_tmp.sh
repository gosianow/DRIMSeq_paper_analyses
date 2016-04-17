#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout


##############################################################################
# Run
##############################################################################

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout


workers=1

for filter_method in 'filter0' 
do
for simulation in 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 
  do
    
    echo "${simulation}_${count_method}_${filter_method}"
    
    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common' disp_prior_df=0.1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common1.Rout  
    
    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='trended' disp_prior_df=1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_trended1.Rout  

  done
done
done



RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout


workers=1

for filter_method in 'filter0'
do
for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common' disp_prior_df=0.1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common2.Rout  
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='trended' disp_prior_df=1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_trended2.Rout  

  done
done
done


RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout


workers=1

for filter_method in 'filter1'
do
for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common' disp_prior_df=0.1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common3.Rout  
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='trended' disp_prior_df=1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_trended3.Rout  

  done
done
done



RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout


workers=1

for filter_method in 'filter1'
do
for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common' disp_prior_df=0.1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common4.Rout  
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='trended' disp_prior_df=1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_trended4.Rout  

  done
done
done



RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout


workers=1

for filter_method in 'filter2'
do
for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common' disp_prior_df=0.1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common5.Rout  
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='trended' disp_prior_df=1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_trended5.Rout  

  done
done
done



RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout


workers=1

for filter_method in 'filter3'
do
for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common' disp_prior_df=0.1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common6.Rout  
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='trended' disp_prior_df=1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_trended6.Rout  

  done
done
done












