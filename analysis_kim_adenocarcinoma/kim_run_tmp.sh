#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Human/Ensembl_GRCh37.71


###############################################################################
### DRIMSeq
###############################################################################


RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Human/Ensembl_GRCh37.71


for model in 'model_full'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"


    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common' disp_prior_df=0.1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_common1.Rout
    
    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='trended' disp_prior_df=1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_trended1.Rout

  done
done



RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Human/Ensembl_GRCh37.71


for model in 'model_null_normal1'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"


    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common' disp_prior_df=0.1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_common2.Rout
    
    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='trended' disp_prior_df=1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_trended2.Rout

  done
done



RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Human/Ensembl_GRCh37.71


for model in 'model_null_normal2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"


    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common' disp_prior_df=0.1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_common3.Rout
    
    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='trended' disp_prior_df=1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_trended3.Rout

  done
done



RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Human/Ensembl_GRCh37.71


for model in 'model_null_tumor1'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"


    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common' disp_prior_df=0.1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_common4.Rout
    
    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='trended' disp_prior_df=1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_trended4.Rout

  done
done



RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Human/Ensembl_GRCh37.71


for model in 'model_null_tumor2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"


    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common' disp_prior_df=0.1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_common5.Rout
    
    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='trended' disp_prior_df=1" $RCODE/kim_drimseq_0_3_3_run.R $ROUT/kim_drimseq_0_3_3_run_${model}_${count_method}_grid_trended5.Rout

  done
done









