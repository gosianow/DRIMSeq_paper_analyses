#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters

# mkdir $ROUT

## Run R scripts



n=3
nm="c(1000,10000)"
param_pi_path="c('$DMPARAMS/kim_kallisto/prop_q15_kim_kallisto_overall.txt','$DMPARAMS/kim_kallisto/prop_q20_kim_kallisto_fcutoff.txt')"
out_name_plots='all_'

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' n=${n} nm=${nm} nd=0 param_pi_path=${param_pi_path} param_gamma_path='$DMPARAMS/kim_kallisto/disp_common_kim_kallisto.txt' out_name_plots='${out_name_plots}'" $RCODE/dispersion_error_filtering_plots_run.R $ROUT/dispersion_error_filtering_plots_run.Rout
     
    

#############################################################################
### Test
#############################################################################



#############################################################################
### Individual run
#############################################################################






















