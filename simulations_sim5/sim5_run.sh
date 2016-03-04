#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout

mkdir $ROUT

##############################################################################
# Run
##############################################################################


workers=5

for filter_method in 'filter0' 'filter1' 'filter2' 'filter3'
do
for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=TRUE results_common=TRUE disp_mode='grid' disp_moderation='none'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_none.Rout
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common.Rout  
    

  done
done
done



############################
# Individual runs
############################



##############################################################################
### Histograms of features
##############################################################################


for filter_method in 'filter0' 'filter1' 'filter2' 'filter3'
do
for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}'" $RCODE/sim5_histograms_of_features.R $ROUT/sim5_histograms_of_features.Rout

  done
done
done


##############################################################################
### Comparison 
##############################################################################

for filter_method in 'filter0' 'filter1' 'filter2' 'filter3'
do

for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}'" $RCODE/sim5_drimseq_comparison_run.R $ROUT/sim5_drimseq_comparison_run.Rout

  done
done
done


##############################################################################
# Combined plots
##############################################################################

for filter_method in 'filter0' 'filter1' 'filter2' 'filter3'
do
  
  echo "${filter_method}"
# All

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull') count_method_list=c('kallisto','kallistofiltered5','kallistoprefiltered5','htseq','htseqprefiltered5') filter_method='${filter_method}' name='_all' legend_nrow=1 pdf_width=18 pdf_height=10" $RCODE/sim5_drimseq_comparison_combined.R $ROUT/sim5_drimseq_comparison_combined.Rout


# Main

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull') count_method_list=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5') filter_method='${filter_method}' name='' legend_nrow=1 pdf_width=14 pdf_height=8" $RCODE/sim5_drimseq_comparison_combined.R $ROUT/sim5_drimseq_comparison_combined.Rout

done



##############################################################################
# Plots with different filtering
##############################################################################



R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull') count_method_list=c('kallisto','htseq') filter_method_list=c('filter0','filter1','filter2','filter3') prefilter_method_list=NULL name='' legend_nrow=2 pdf_width=10 pdf_height=8" $RCODE/sim5_drimseq_comparison_filtering.R $ROUT/sim5_drimseq_comparison_filtering.Rout




R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull') count_method_list=c('kallisto','htseq') filter_method_list=c('filter0','filter1','filter2','filter3') prefilter_method_list=c('kallistoprefiltered5','htseqprefiltered5') name='_prefilt' legend_nrow=2 pdf_width=10 pdf_height=8" $RCODE/sim5_drimseq_comparison_filtering.R $ROUT/sim5_drimseq_comparison_filtering.Rout





###################################################################################################################

### Using the F-test

###################################################################################################################

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout

mkdir $ROUT


##############################################################################
# Run using the F test - R32dev
##############################################################################


workers=5

for filter_method in 'filter0'
do
for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"

    R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=TRUE results_common=FALSE disp_mode='grid' disp_moderation='none'" $RCODE/sim5_drimseq_f_run.R $ROUT/sim5_drimseq_f_run_${simulation}_${count_method}_${filter_method}_grid_none.Rout
    
    R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common'" $RCODE/sim5_drimseq_f_run.R $ROUT/sim5_drimseq_f_run_${simulation}_${count_method}_${filter_method}_grid_common.Rout  
    

  done
done
done




##############################################################################
### Comparison 
##############################################################################


for filter_method in 'filter0'
do

for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${simulation}_${count_method}_${filter_method}"

    R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}'" $RCODE/sim5_drimseq_f_comparison_run.R $ROUT/sim5_drimseq_f_comparison_run.Rout

  done
done
done



##############################################################################
# Combined plots
##############################################################################

for filter_method in 'filter0'
do
  
  echo "${filter_method}"

# Main

R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull') count_method_list=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5') filter_method='${filter_method}' name='' legend_nrow=1 pdf_width=14 pdf_height=8" $RCODE/sim5_drimseq_f_comparison_combined.R $ROUT/sim5_drimseq_f_comparison_combined.Rout

done




