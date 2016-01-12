## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout

mkdir $ROUT


### Run R scripts

# model full

for model in 'model_full'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=TRUE results_common=TRUE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_${model}_${count_method}_grid_none.Rout
    
    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_${model}_${count_method}_grid_common.Rout
    
    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list=c('optimize','optim','constrOptim') disp_moderation_list=c('none','none','none')" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_${model}_${count_method}_optim.Rout

  done
done


# other models

for model in 'model_null_normal1' 'model_null_normal2' 'model_null_tumor1' 'model_null_tumor2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_${model}_${count_method}_grid_none.Rout
    
    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_${model}_${count_method}_grid_none.Rout

  done
done


###############################################################################
### Individual runs
###############################################################################



for model in 'model_full' 'model_null_normal1' 'model_null_normal2' 'model_null_tumor1' 'model_null_tumor2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_add_d.Rout

  done
done
























