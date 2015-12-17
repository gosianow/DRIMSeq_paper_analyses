## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout

mkdir $ROUT

##############################################################################
# Run
##############################################################################

filter_method='filter1'
workers=5

for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
  do
    
    echo "${simulation}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=TRUE results_common=TRUE disp_mode='grid' disp_moderation='none'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_none.Rout
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common.Rout  
    
    # R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='optimize' disp_moderation='none'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_optimize.Rout
    
    # R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='optim' disp_moderation='none'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_optim.Rout
    
    # R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='constrOptim' disp_moderation='none'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_constrOptim.Rout

  done
done



filter_method='filter3'
workers=20

for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
  do
    
    echo "${simulation}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=TRUE results_common=TRUE disp_mode='grid' disp_moderation='none'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_none.Rout
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common.Rout  
    
  done
done



############################
# Individual runs
############################



##############################################################################
### Histograms of features
##############################################################################

filter_method='filter0'

for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5'
  do
    
    echo "${simulation}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}'" $RCODE/sim5_histograms_of_features.R $ROUT/sim5_histograms_of_features.Rout

  done
done


##############################################################################
### Comparison 
##############################################################################

filter_method='filter1'

for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5'
  do
    
    echo "${simulation}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}'" $RCODE/sim5_drimseq_comparison_run.R $ROUT/sim5_drimseq_comparison_run.Rout

  done
done



##############################################################################
# Combined plots
##############################################################################

filter_method='filter0'
# All

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull') count_method_list=c('kallisto','kallistofiltered5','kallistoprefiltered5','htseq','htseqprefiltered5') filter_method='${filter_method}' name='_all' legend_nrow=1 pdf_width=18 pdf_height=10" $RCODE/sim5_drimseq_comparison_combined.R $ROUT/sim5_drimseq_comparison_combined.Rout


# Main

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull') count_method_list=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5') filter_method='${filter_method}' name='' legend_nrow=1 pdf_width=14 pdf_height=8" $RCODE/sim5_drimseq_comparison_combined.R $ROUT/sim5_drimseq_comparison_combined.Rout
























