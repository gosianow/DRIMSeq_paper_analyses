## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout

# mkdir $ROUT

## Run R scripts

# null_normal2

R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='kallisto' model='model_null_normal2'" $RCODE/kim_dexseq_run.R $ROUT/kim_dexseq_run_null_normal2_kallisto.Rout

R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='htseq' model='model_null_normal2'" $RCODE/kim_dexseq_run.R $ROUT/kim_dexseq_run_null_normal2_htseq.Rout


# null_tumor2

R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='kallisto' model='model_null_tumor2'" $RCODE/kim_dexseq_run.R $ROUT/kim_dexseq_run_null_tumor2_kallisto.Rout

R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='htseq' model='model_null_tumor2'" $RCODE/kim_dexseq_run.R $ROUT/kim_dexseq_run_null_tumor2_htseq.Rout










