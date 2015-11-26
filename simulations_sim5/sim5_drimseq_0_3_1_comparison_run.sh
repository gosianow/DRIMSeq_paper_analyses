## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout

# mkdir $ROUT

### Run R scripts

### Comparison plots

for i in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for j in 'kallisto' 'htseq' 'htseqprefiltered15' 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5'
  do
    
    echo "${i}_${j}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${j}' simulation='${i}'" $RCODE/sim5_drimseq_0_3_1_comparison_run.R $ROUT/sim5_drimseq_0_3_1_comparison_run.Rout

  done
done


### Histograms of features

for i in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
do 
  for j in 'kallisto' 'htseq' 'htseqprefiltered15' 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5'
  do
    
    echo "${i}_${j}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${j}' simulation='${i}'" $RCODE/sim5_histograms_of_features.R $ROUT/sim5_histograms_of_features.Rout

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





























