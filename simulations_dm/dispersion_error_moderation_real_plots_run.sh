#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters

# mkdir $ROUT

## Run R scripts


for n in 3
do
  
for data_name in 'kim' 'brooks'
do 
  
for count_method in 'kallisto' 'htseq'
do
    

      echo "${data_name}_${count_method}_n${n}"
      
    
    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/dispersion_error_moderation_real_plots_run.R $ROUT/dispersion_error_moderation_real_plots_run_${data_name}_${count_method}_n${n}.Rout
    

done
done
done




##############################################################################
### Test
##############################################################################



##############################################################################
### Individual runs
##############################################################################

