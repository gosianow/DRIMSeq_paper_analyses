## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout

# mkdir $ROUT

## Run R scripts

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/dispersion_error_combined_plots_run.R $ROUT/dispersion_error_combined_plots_run.Rout







