#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_sim5
RWD=/home/gosia/multinomial_project/simulations_sim5
ROUT=$RWD/Rout

mkdir -p $ROUT


##############################
### Colors
##############################

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD'" $RCODE/colors.R $ROUT/colors.Rout

###############################################################################
### Run only for hsapiens_withde_nonull because the filtering analysis were missing in Simulation5_Charlotte
###############################################################################

# I have coppied from /home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/hsapiens/reference_files/ into hsapiens_withde_nonull/0_annotation/ the following:
# TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa, Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf, KallistoIndex/TranscriptID_conversion.txt
# And /home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/hsapiens/with_diffexpression/non_null_simulation/1_reads

#######################################
# Kallisto
#######################################

### Obtain kallisto transcrit abundance based on the original transcript catalog
R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD/hsapiens_withde_nonull' cDNA_fasta='0_annotation/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa' gtf_path='0_annotation/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' conversion_path='0_annotation/KallistoIndex/TranscriptID_conversion.txt' out_dir='2_counts/kallisto'" $RCODE/sim5_kallisto.R $ROUT/sim5_kallisto.Rout
tail $ROUT/sim5_kallisto.Rout

### Filter out the lowly expressed transcripts (from GTF, FA and kallisto expression)
R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD/hsapiens_withde_nonull' fa_path='0_annotation/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa' gtf_path='0_annotation/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' conversion_path='0_annotation/KallistoIndex/TranscriptID_conversion.txt' kallisto_dir='2_counts/kallisto' out_dir='2_counts/kallisto_txfilt_5'" $RCODE/sim5_kallisto_filter_gtf_fa.R $ROUT/sim5_kallisto_filter_gtf_fa.Rout
tail $ROUT/sim5_kallisto_filter_gtf_fa.Rout

### Obtain kallisto transcrit abundance based on the prefiltered transcript catalog
R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD/hsapiens_withde_nonull' cDNA_fasta='0_annotation/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.dna.primary_assembly_kallistoest_atleast5.fa' gtf_path='0_annotation/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_kallistoest_atleast5.gtf' conversion_path='0_annotation/KallistoIndex/TranscriptID_conversion.txt' out_dir='2_counts/INCOMPLETE_KALLISTOEST/kallisto_kallistoest_atleast5'" $RCODE/sim5_kallisto.R $ROUT/sim5_kallisto.Rout
tail $ROUT/sim5_kallisto.Rout

#######################################
# HTSeq counts
#######################################

### counts for the original transcript catalog
R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD/hsapiens_withde_nonull' gtf_path='0_annotation/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf out_dir=2_counts/dexseq_nomerge" $RCODE/sim5_htseq.R $ROUT/sim5_htseq.Rout
tail $ROUT/sim5_htseq.Rout

### counts for the prefiltered catalog
R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD/hsapiens_withde_nonull' gtf_path='0_annotation/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_kallistoest_atleast5.gtf'  out_dir='2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast5'" $RCODE/sim5_htseq.R $ROUT/sim5_htseq.Rout
tail $ROUT/sim5_htseq.Rout

#######################################
# Run DEXSeq
#######################################

for count_method in 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
do

  R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD/hsapiens_withde_nonull' workers=5 count_method='${count_method}'" $RCODE/sim5_dexseq_run.R $ROUT/sim5_dexseq_run.Rout
  tail $ROUT/sim5_dexseq_run.Rout

done


##############################################################################
# Run DRIMSeq
##############################################################################


workers=1

for filter_method in 'filter0' 'filter1' 'filter2' 'filter3'
do
  for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
  do
    for count_method in 'kallisto' 'htseq' 'kallistoprefiltered5' 'htseqprefiltered5' 'kallistofiltered5'
    do

      echo "${simulation}_${count_method}_${filter_method}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=TRUE results_common=TRUE disp_mode='grid' disp_moderation='none' disp_prior_df=0.1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_none.Rout

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='common' disp_prior_df=0.1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_common.Rout

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' dispersion_common=FALSE results_common=FALSE disp_mode='grid' disp_moderation='trended' disp_prior_df=1" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run_${simulation}_${count_method}_${filter_method}_grid_trended.Rout

    done
  done
done



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

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison'" $RCODE/sim5_histograms_of_features.R $ROUT/sim5_histograms_of_features.Rout

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
    for count_method in 'kallistofiltered5' 'htseqprefiltered5' 'kallistoprefiltered5' 'kallisto' 'htseq'
    do

      echo "${simulation}_${count_method}_${filter_method}"

      R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R' method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison'" $RCODE/sim5_drimseq_comparison.R $ROUT/sim5_drimseq_comparison_run.Rout

      tail $RCODE/sim5_drimseq_comparison.R $ROUT/sim5_drimseq_comparison_run.Rout

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

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull') count_method_list=c('kallisto','kallistofiltered5','kallistoprefiltered5','htseq','htseqprefiltered5') filter_method='${filter_method}' CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R' name='_all' legend_nrow=1 pdf_width=18 pdf_height=10 method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison' strip_text_size=12 text_size=16 legend_size=18" $RCODE/sim5_drimseq_comparison_combined.R $ROUT/sim5_drimseq_comparison_combined.Rout

  tail -v $ROUT/sim5_drimseq_comparison_combined.Rout


  # Main

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull') count_method_list=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5') filter_method='${filter_method}' CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R' name='' legend_nrow=1 pdf_width=14 pdf_height=8 method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison' strip_text_size=12 text_size=16 legend_size=15" $RCODE/sim5_drimseq_comparison_combined.R $ROUT/sim5_drimseq_comparison_combined.Rout

  tail -v $ROUT/sim5_drimseq_comparison_combined.Rout

done




##############################################################################
# Plots with different filtering
##############################################################################



R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull') count_method_list=c('kallisto','htseq') filter_method_list=c('filter0','filter1','filter2','filter3') prefilter_method_list=NULL name='' legend_nrow=3 pdf_width=12 pdf_height=9 method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison' strip_text_size=14 text_size=16 legend_size=14" $RCODE/sim5_drimseq_comparison_filtering.R $ROUT/sim5_drimseq_comparison_filtering.Rout


### Prefiltering plotted as a separate method

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull') count_method_list=c('kallisto','htseq') filter_method_list=c('filter0','filter1','filter2','filter3') prefilter_method_list=c('kallistoprefiltered5','htseqprefiltered5') name='_prefilt' legend_nrow=3 pdf_width=12 pdf_height=9 method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison' strip_text_size=14 text_size=16 legend_size=14" $RCODE/sim5_drimseq_comparison_filtering.R $ROUT/sim5_drimseq_comparison_filtering.Rout




##############################################################################
# Run auto moderation diagnostics
##############################################################################



workers=10

for filter_method in 'filter0'
do
  for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull' 'hsapiens_withde_nonull'
  do
    for count_method in 'kallistofiltered5' 'htseqprefiltered5' 'kallisto' 'htseq' 'kallistoprefiltered5'
    do

      echo "${simulation}_${count_method}_${filter_method}"

      R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=${workers} count_method='${count_method}' simulation='${simulation}' filter_method='${filter_method}' method_out='drimseq_0_3_3' dmDS_auto_moderation_diagnostics_function_path='/home/gosia/R/drimseq_paper/help_functions/dmDS_auto_moderation_diagnostics.R'" $RCODE/sim5_drimseq_auto_moderation.R $ROUT/sim5_drimseq_auto_moderation_${simulation}_${count_method}_${filter_method}.Rout

      tail $ROUT/sim5_drimseq_auto_moderation_${simulation}_${count_method}_${filter_method}.Rout

    done
  done
done
