#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq
ROUT=$RWD/Rout

mkdir $ROUT

## Run R scripts

### Kim 
R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' main_data_dir='/home/Shared/data/seq/kim_adenocarcinoma' data_name='kim' count_method='kallisto'" $RCODE/dm_parameters_run.R $ROUT/dm_parameters_run_kim_kallisto.Rout


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' main_data_dir='/home/Shared/data/seq/kim_adenocarcinoma' data_name='kim' count_method='htseq'" $RCODE/dm_parameters_run.R $ROUT/dm_parameters_run_kim_htseq.Rout

### Brooks
R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' main_data_dir='/home/Shared/data/seq/brooks_pasilla' data_name='brooks' count_method='kallisto'" $RCODE/dm_parameters_run.R $ROUT/dm_parameters_run_brooks_kallisto.Rout


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' main_data_dir='/home/Shared/data/seq/brooks_pasilla' data_name='brooks' count_method='htseq'" $RCODE/dm_parameters_run.R $ROUT/dm_parameters_run_brooks_htseq.Rout



