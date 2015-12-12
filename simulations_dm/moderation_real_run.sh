#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters

# mkdir $ROUT


##############################################################################
### Run
##############################################################################


workers=5

for n in 3 6
do
  
for data_name in 'kim' 'brooks'
do 
  
for count_method in 'kallisto' 'htseq'
do
    
    for run in {1..20}
    do
      
      echo "${data_name}_${count_method}_n${n}_${run}"
      
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=5000 disp_prior_df=c(0,c(1,3,5)/10000,c(1,3,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,3,5,10) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/moderation_real_run.R $ROUT/moderation_real_run_${data_name}_${count_method}_n${n}.Rout
    
  done

done
done
done



################################
### Test
################################

for n in 3
do
  
for data_name in 'kim'
do 
  
for count_method in 'kallisto'
do
    
    for run in {1..3}
    do
      
      echo "${data_name}_${count_method}_n${n}_${run}"
      
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=4 sim_name='test_' run='run${run}' m=100 disp_prior_df=c(0,c(1)/100,c(1)/10,1,10) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/moderation_real_run.R $ROUT/moderation_real_run_${data_name}_${count_method}_n${n}.Rout
    
  done

done
done
done


################################
### Individual runs
################################

workers=5

for n in 3 6
do
  
for data_name in 'brooks'
do 
  
for count_method in 'kallisto' 'htseq'
do
    
    for run in {1..20}
    do
      
      echo "${data_name}_${count_method}_n${n}_${run}"
      
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=5000 disp_prior_df=c(0,c(1,3,5)/10000,c(1,3,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,3,5,10) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/moderation_real_run.R $ROUT/moderation_real_run_${data_name}_${count_method}_n${n}.Rout
    
  done

done
done
done


workers=5

for n in 3 6
do
  
for data_name in 'kim'
do 
  
for count_method in 'kallisto' 'htseq'
do
    
    for run in {1..20}
    do
      
      echo "${data_name}_${count_method}_n${n}_${run}"
      
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=5000 disp_prior_df=c(0,c(1,3,5)/10000,c(1,3,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,3,5,10) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/moderation_real_run.R $ROUT/moderation_real_run_${data_name}_${count_method}_n${n}.Rout
    
  done

done
done
done


##############################################################################
### Plots
##############################################################################



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




