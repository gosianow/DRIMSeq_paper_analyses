## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters

mkdir $ROUT

## Run R scripts

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='s4_' r=1 m=500 n=6 nm=1000 disp_prior_df=c(seq(0,1,by=0.1),2,3) param_pi_path='$DMPARAMS/prop_q3_uniform.txt' param_gamma_path='$DMPARAMS/disp_genewise_kim_kallisto_lognormal.txt'" $RCODE/dispersion_error_moderation_run.R $ROUT/dispersion_error_moderation_run_test.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='s4_' r=1 m=500 n=6 nm=1000 disp_prior_df=c(seq(0,1,by=0.1),2,3) param_pi_path='$DMPARAMS/prop_q10_uniform.txt' param_gamma_path='$DMPARAMS/disp_genewise_kim_kallisto_lognormal.txt'" $RCODE/dispersion_error_moderation_run.R $ROUT/dispersion_error_moderation_run_test.Rout
