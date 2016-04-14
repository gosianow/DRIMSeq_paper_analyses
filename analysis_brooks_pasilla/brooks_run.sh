#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_brooks_pasilla
RWD=/home/Shared/data/seq/brooks_pasilla
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Drosophila/Ensembl70

mkdir $ROUT

## Run R scripts


###############################################################################
### Download data + Mapping with tophat2
###############################################################################

### ! First save SraRunInfo.csv file in "$RWD/3_metadata" directory !


R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' DNA_fasta='$ANNOTATION/genome/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa'" $RCODE/brooks_download.R $ROUT/brooks_download.Rout



###############################################################################
### Kallisto
###############################################################################


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf_path='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' cDNA_fasta='$ANNOTATION/cDNA/Drosophila_melanogaster.BDGP5.70.cdna.all.fa'" $RCODE/brooks_kallisto.R $ROUT/brooks_kallisto.Rout



R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf_path='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf'" $RCODE/brooks_kallisto_filter_gtf.R $ROUT/brooks_kallisto_filter_gtf.Rout


###############################################################################
### DEXSeq
###############################################################################


### Index BAM files

for i in 'GSM461176' 'GSM461177' 'GSM461178' 'GSM461179' 'GSM461180' 'GSM461181' 'GSM461182'
 do 
  samtools index $RWD/1_reads/tophat_2.0.14/${i}/accepted_hits.bam
done


### Run HTSeq

R214 CMD BATCH --no-save --no-restore "--args gtf='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' count_method='htseq'" $RCODE/brooks_htseq.R $ROUT/brooks_htseq_htseq.Rout


R214 CMD BATCH --no-save --no-restore "--args gtf='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70_kallistoest_atleast5.gtf' count_method='htseqprefiltered5'" $RCODE/brooks_htseq.R $ROUT/brooks_htseq_htseqprefiltered5.Rout


### Run DEXSeq

for model in 'model_full' 'model_full_glm' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=10 count_method='${count_method}' model='${model}'" $RCODE/brooks_dexseq_run.R $ROUT/brooks_dexseq_run_${model}_${count_method}.Rout

  done
done



###############################################################################
### DRIMSeq
###############################################################################


for model in 'model_full' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none' disp_prior_df=0.1" $RCODE/brooks_drimseq_0_3_3_run.R $ROUT/brooks_drimseq_0_3_3_run_${model}_${count_method}_grid_none.Rout
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common' disp_prior_df=0.1" $RCODE/brooks_drimseq_0_3_3_run.R $ROUT/brooks_drimseq_0_3_3_run_${model}_${count_method}_grid_common.Rout
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='trended' disp_prior_df=1" $RCODE/brooks_drimseq_0_3_3_run.R $ROUT/brooks_drimseq_0_3_3_run_${model}_${count_method}_grid_trended.Rout

  done
done


##############################
### Colors
##############################

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/drimseq_0_3_3_comparison'" $RCODE/colors.R $ROUT/colors.Rout

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/drimseq_0_3_3_positive_controls'" $RCODE/colors.R $ROUT/colors.Rout


##############################
### DRIMSeq comparison
##############################

### Plot venn diagrams, upset plots; Update dispersion plots; Create a summary file with numbers of positive genes

for model in 'model_full' 'model_full_glm' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${model}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' model='${model}' method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison'" $RCODE/brooks_drimseq_0_3_3_comparison.R $ROUT/brooks_drimseq_0_3_3_comparison_run.Rout

  done
done


### Barplots of the number of all and DS genes + overlaps

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' comparison_out='drimseq_0_3_3_comparison' keep_methods=c('dexseq','drimseq_genewise_grid_none','drimseq_genewise_grid_common','drimseq_genewise_grid_trended')" $RCODE/brooks_drimseq_0_3_3_comparison_summary.R $ROUT/brooks_drimseq_0_3_3_summary.Rout


### Plots of the overlap versus number of top ranked genes + CAT plots

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD'  count_methods=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5') models=c('model_full','model_full_paired') method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison' Overlaps_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotOverlaps.R' CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R'" $RCODE/brooks_drimseq_0_3_3_comparison_plots.R $ROUT/brooks_drimseq_0_3_3_comparison_plots.Rout


### Plot overlaps between models

for ds_method in 'drimseq_genewise_grid_none' 'drimseq_genewise_grid_common' 'drimseq_genewise_grid_trended'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${ds_method}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' ds_method='${ds_method}' model_list=c('model_full','model_full_paired') method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison'" $RCODE/brooks_drimseq_0_3_3_comparison_models.R $ROUT/brooks_drimseq_0_3_3_comparison_models.Rout

  done
done


for ds_method in 'dexseq'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${ds_method}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' ds_method='${ds_method}' model_list=c('model_full','model_full_paired','model_full_glm') method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison'" $RCODE/brooks_drimseq_0_3_3_comparison_models.R $ROUT/brooks_drimseq_0_3_3_comparison_models.Rout

  done
done


##############################
### DRIMSeq positive controls
##############################

### ! Save brooks_validated_genes.txt in "$RWD/5_validation" directory !

### Add Ensemble gene IDs to the brooks_validated_genes.txt file

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf_path='/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf'" $RCODE/brooks_validated_genes.R $ROUT/brooks_validated_genes.Rout


### Make DRIMSeq plots of proportions for the validated genes; Create summary files

for model in 'model_full' 'model_full_paired' 'model_full_glm'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${model}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' valid_path='5_validation/brooks_validated_genes.txt' count_method='${count_method}' model='${model}' method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_positive_controls'" $RCODE/brooks_drimseq_0_3_3_positive_controls.R $ROUT/brooks_drimseq_0_3_3_positive_controls_run.Rout

  done
done


### Plot tables with adj p-values for validated genes - validation_summary.pdf; Plot coverage and annotations for the validated genes

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_gtf='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' path_gtf_filtered='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70_kallistoest_atleast5.gtf' valid_path='5_validation/brooks_validated_genes.txt' method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_positive_controls' keep_methods=c('dexseq','drimseq_genewise_grid_none','drimseq_genewise_grid_common','drimseq_genewise_grid_trended')" $RCODE/brooks_drimseq_0_3_3_positive_controls_summary.R $ROUT/brooks_drimseq_0_3_3_positive_controls_summary.Rout



### Plots of the number of validated genes versus number of genes detected as DS


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' valid_path='5_validation/brooks_validated_genes.txt' count_methods=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5') models=c('model_full','model_full_paired','model_full_glm') method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_positive_controls' ROC_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotROCx.R'" $RCODE/brooks_drimseq_0_3_3_positive_controls_plots.R $ROUT/brooks_drimseq_0_3_3_positive_controls_plots.Rout





















