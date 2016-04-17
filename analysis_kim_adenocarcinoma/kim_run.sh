#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Human/Ensembl_GRCh37.71

mkdir $ROUT

## Run R scripts

###############################################################################
### Download FASTQ files and bam files from insilicodb
###############################################################################

### ! First save SraRunInfo.csv file in "$RWD/3_metadata" directory !
### ! First save Malgorzata Nowicka2014-11-04GSE37764.csv file in "$RWD/3_metadata" directory !


R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/kim_download.R $ROUT/kim_download.Rout


###############################################################################
### Kallisto
###############################################################################


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf_path='$ANNOTATION/gtf/Homo_sapiens.GRCh37.71.gtf' cDNA_fasta='$ANNOTATION/cDNA/Homo_sapiens.GRCh37.71.cdna.all.fa'" $RCODE/kim_kallisto.R $ROUT/kim_kallisto.Rout



R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf_path='$ANNOTATION/gtf/Homo_sapiens.GRCh37.71.gtf'" $RCODE/kim_kallisto_filter_gtf.R $ROUT/kim_kallisto_filter_gtf.Rout



###############################################################################
### DEXSeq
###############################################################################


### Index BAM files

for i in 'GSM927308' 'GSM927310' 'GSM927312' 'GSM927314' 'GSM927316' 'GSM927318' 'GSM927309' 'GSM927311' 'GSM927313' 'GSM927315' 'GSM927317' 'GSM927319'
 do 
  samtools index $RWD/1_reads/tophat_insilicodb/${i}/accepted_hits.bam
done


### Run HTSeq


R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf='$ANNOTATION/gtf/Homo_sapiens.GRCh37.71.gtf' count_method='htseq'" $RCODE/kim_htseq.R $ROUT/kim_htseq_htseq.Rout


R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf='$ANNOTATION/gtf/Homo_sapiens.GRCh37.71_kallistoest_atleast5.gtf' count_method='htseqprefiltered5'" $RCODE/kim_htseq.R $ROUT/kim_htseq_htseqprefiltered5.Rout




### Run DEXSeq

for model in 'model_full' 'model_full_glm' 'model_null_normal1' 'model_null_normal2' 'model_null_tumor1' 'model_null_tumor2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do

    echo "${model}_${count_method}"

    R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}'" $RCODE/kim_dexseq_run.R $ROUT/kim_dexseq_run_${model}_${count_method}.Rout

  done
done



###############################################################################
### DRIMSeq
###############################################################################


for model in 'model_full' 'model_null_normal1' 'model_null_normal2' 'model_null_tumor1' 'model_null_tumor2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none' disp_prior_df=0.1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_none.Rout
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common' disp_prior_df=0.1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_common.Rout
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='trended' disp_prior_df=1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_trended.Rout

  done
done



##############################
### Colors
##############################

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/drimseq_0_3_3_comparison'" $RCODE/colors.R $ROUT/colors.Rout


##############################
### DRIMSeq comparison
##############################

### Plot venn diagrams, upset plots; Update dispersion plots; Create a summary file with numbers of positive genes

for model in 'model_full' 'model_full_glm' 'model_null_normal1' 'model_null_normal2' 'model_null_tumor1' 'model_null_tumor2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${model}_${count_method}"

    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' model='${model}' method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison'" $RCODE/kim_drimseq_0_3_3_comparison.R $ROUT/kim_drimseq_0_3_3_comparison_run.Rout

  done
done



### Barplots of the number of all and DS genes + overlaps

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' comparison_out='drimseq_0_3_3_comparison' keep_methods=c('dexseq','drimseq_genewise_grid_none','drimseq_genewise_grid_common','drimseq_genewise_grid_trended')" $RCODE/kim_drimseq_0_3_3_comparison_summary.R $ROUT/kim_drimseq_0_3_3_comparison_summary.Rout



### Plot the overlap versus number of top ranked genes + CAT plots

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD'  count_methods=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5') models=c('model_full') method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison' Overlaps_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotOverlaps.R' CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R'" $RCODE/kim_drimseq_0_3_3_comparison_plots.R $ROUT/kim_drimseq_0_3_3_comparison_plots.Rout


### Plot overlaps between models

for ds_method in 'dexseq'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${ds_method}_${count_method}"

    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' ds_method='${ds_method}' model_list=c('model_full','model_full_glm') method_out='drimseq_0_3_3' comparison_out='drimseq_0_3_3_comparison'" $RCODE/kim_drimseq_0_3_3_comparison_models.R $ROUT/kim_drimseq_0_3_3_comparison_models.Rout

  done
done












