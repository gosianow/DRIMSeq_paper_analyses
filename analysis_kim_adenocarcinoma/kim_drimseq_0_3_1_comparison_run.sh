## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout

mkdir $ROUT

## Run R scripts

# full

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='kallisto' model='model_full'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='htseq' model='model_full'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

# null 1

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='kallisto' model='model_null_normal1'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='htseq' model='model_null_normal1'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='kallisto' model='model_null_tumor1'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='htseq' model='model_null_tumor1'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

# null 2

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='kallisto' model='model_null_normal2'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='htseq' model='model_null_normal2'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='kallisto' model='model_null_tumor2'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='htseq' model='model_null_tumor2'" $RCODE/kim_drimseq_0_3_1_comparison_run.R $ROUT/kim_drimseq_0_3_1_comparison_run.Rout







