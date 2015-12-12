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

############################
# Individual runs
############################

filter_method='filter1'
workers=5

for simulation in 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
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



##############################################################################
### Histograms of features
##############################################################################

for i in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for j in 'kallisto' 'htseq' 'htseqprefiltered15' 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5'
  do
    
    echo "${i}_${j}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${j}' simulation='${i}'" $RCODE/sim5_histograms_of_features.R $ROUT/sim5_histograms_of_features.Rout

  done
done


##############################################################################
### Comparison 
##############################################################################


for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'htseqprefiltered15' 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5'
  do
    
    echo "${simulation}_${count_method}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}'" $RCODE/sim5_drimseq_comparison_run.R $ROUT/sim5_drimseq_comparison_run.Rout

  done
done



##############################################################################
# Combined plots
##############################################################################

# All

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull') count_method_list=c('kallisto','kallistofiltered5','kallistoprefiltered5','htseq','htseqprefiltered5','htseqprefiltered15') name='_all' legend_nrow=1 pdf_width=18 pdf_height=10" $RCODE/sim5_drimseq_0_3_1_comparison_combined.R \
$ROUT/sim5_drimseq_0_3_1_comparison_combined.Rout


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull') count_method_list=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5') name='' legend_nrow=1 pdf_width=14 pdf_height=8" $RCODE/sim5_drimseq_0_3_1_comparison_combined.R \
$ROUT/sim5_drimseq_0_3_1_comparison_combined.Rout




