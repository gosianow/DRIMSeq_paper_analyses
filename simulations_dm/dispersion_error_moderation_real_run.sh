## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters

# mkdir $ROUT

## Run R scripts

# lognormal dispersion from kallisto: disp_genewise_kim_kallisto_lognormal.txt






##############################################################################
### Test
##############################################################################



for n in 3
do
  
  echo "n${n}"
  
  R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=10 sim_name='test3_' r=5 m=10000 disp_prior_df=c(c(1,3,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,3,5,10) param_pi_path='$DMPARAMS/prop_kim_kallisto.txt' param_gamma_path='$DMPARAMS/disp_genewise_kim_kallisto_lognormal.txt' param_nm_path='$DMPARAMS/nm_kim_kallisto_lognormal.txt' param_nd_path='$DMPARAMS/nd_common_kim_kallisto.txt'" $RCODE/dispersion_error_moderation_real_run.R $ROUT/dispersion_error_moderation_real_run_n${n}.Rout


done


##############################################################################
### Individual runs
##############################################################################











