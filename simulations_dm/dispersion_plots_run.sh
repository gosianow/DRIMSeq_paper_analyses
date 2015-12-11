#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout

# mkdir $ROUT

## Run R scripts

### Combined plots of error

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/dispersion_error_combined_plots_run.R $ROUT/dispersion_error_combined_plots_run.Rout

### Combined plots of FP


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/dispersion_fp_combined_plots_run.R $ROUT/dispersion_fp_combined_plots_run.Rout

### Optimization plots


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/dispersion_error_optimization_plots_run.R $ROUT/dispersion_error_optimization_plots_run.Rout



### Moderation plots

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/dispersion_error_moderation_plots_run.R $ROUT/dispersion_error_moderation_plots_run.Rout





### Moderation real plots

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/dispersion_error_moderation_real_plots_run.R $ROUT/dispersion_error_moderation_real_plots_run.Rout










