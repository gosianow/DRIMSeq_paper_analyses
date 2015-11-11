## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout

mkdir $ROUT

## Run R scripts

# full
# kallisto

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='kallisto' model='model_full' dispersion_common=TRUE results_common=TRUE disp_mode_list=c('grid') disp_moderation_list=c('none')" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_full_kallisto_grid_common.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='kallisto' model='model_full' dispersion_common=FALSE results_common=TRUE disp_mode_list=c('optimize','optim','constrOptim') disp_moderation_list=c('none','none','none')" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_full_kallisto.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='kallisto' model='model_full' dispersion_common=FALSE results_common=TRUE disp_mode_list=c('grid') disp_moderation_list=c('common')" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_full_kallisto_grid_common.Rout

# htseq

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='htseq' model='model_full' dispersion_common=TRUE results_common=TRUE disp_mode_list=c('grid','optimize','optim','constrOptim') disp_moderation_list=c('none','none','none','none')" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_full_htseq.Rout


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='htseq' model='model_full' dispersion_common=FALSE results_common=TRUE disp_mode_list=c('grid') disp_moderation_list=c('common')" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_full_htseq_grid_common.Rout



# null_normal1

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=3 count_method='kallisto' model='model_null_normal1' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_null_normal1_kallisto.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=3 count_method='htseq' model='model_null_normal1' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_null_normal1_htseq.Rout


# null_tumor1

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=3 count_method='kallisto' model='model_null_tumor1' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_null_tumor1_kallisto.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=3 count_method='htseq' model='model_null_tumor1' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_null_tumor1_htseq.Rout


# null_normal2

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=3 count_method='kallisto' model='model_null_normal2' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_null_normal1_kallisto.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=3 count_method='htseq' model='model_null_normal2' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_null_normal1_htseq.Rout


# null_tumor2

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=3 count_method='kallisto' model='model_null_tumor2' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_null_tumor1_kallisto.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=3 count_method='htseq' model='model_null_tumor2' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/kim_drimseq_0_3_1_run.R $ROUT/kim_drimseq_0_3_1_run_null_tumor1_htseq.Rout










