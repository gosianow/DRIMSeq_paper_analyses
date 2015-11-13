## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters

# mkdir $ROUT

## Run R scripts

# common dispersion from kallisto & uniform proportions

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=2 m=500 n=3 nm=1000 nd=0 param_pi_path='$DMPARAMS/prop_q3_uniform.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_optimization_run.R $ROUT/dispersion_error_optimization_run1.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=2 m=500 n=3 nm=1000 nd=0 param_pi_path='$DMPARAMS/prop_q10_uniform.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_optimization_run.R $ROUT/dispersion_error_optimization_run2.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=2 m=500 n=3 nm=100 nd=0 param_pi_path='$DMPARAMS/prop_q3_uniform.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_optimization_run.R $ROUT/dispersion_error_optimization_run3.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=2 m=500 n=3 nm=100 nd=0 param_pi_path='$DMPARAMS/prop_q10_uniform.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_optimization_run.R $ROUT/dispersion_error_optimization_run4.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/dispersion_error_optimization_plots_run.R $ROUT/dispersion_error_optimization_plots_run.Rout


# common dispersion from kallisto & descending proportions

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=2 m=500 n=3 nm=1000 nd=0 param_pi_path='$DMPARAMS/prop_q3_kim_kallisto_overall.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_optimization_run.R $ROUT/dispersion_error_optimization_run5.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=2 m=500 n=3 nm=1000 nd=0 param_pi_path='$DMPARAMS/prop_q10_kim_kallisto_overall.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_optimization_run.R $ROUT/dispersion_error_optimization_run6.Rout


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=2 m=500 n=3 nm=100 nd=0 param_pi_path='$DMPARAMS/prop_q3_kim_kallisto_overall.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_optimization_run.R $ROUT/dispersion_error_optimization_run7.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=2 m=500 n=3 nm=100 nd=0 param_pi_path='$DMPARAMS/prop_q10_kim_kallisto_overall.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_optimization_run.R $ROUT/dispersion_error_optimization_run8.Rout

R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/dispersion_error_optimization_plots_run.R $ROUT/dispersion_error_optimization_plots_run.Rout












