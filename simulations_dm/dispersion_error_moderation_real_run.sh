## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters

# mkdir $ROUT

## Run R scripts

# lognormal dispersion from kallisto: disp_genewise_kim_kallisto_lognormal.txt



for n in 3
do
  
for data_name in 'kim'
do 
  
for count_method in 'htseq'
do
    
    for run in {1..20}
    do
      
      echo "${data_name}_${count_method}_n${n}_${run}"
      
    
    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=10 sim_name='run${run}_' m=10000 disp_prior_df=c(0,c(1,3,5)/10000,c(1,3,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,3,5,10) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/dispersion_error_moderation_real_run.R $ROUT/dispersion_error_moderation_real_run_${data_name}_${count_method}_n${n}.Rout
    
  done

done
done
done





for n in 3
do
  
for data_name in 'brooks'
do 
  
for count_method in 'kallisto' 'htseq'
do
  
  for run in {1..20}
  do
    
    echo "${data_name}_${count_method}_n${n}_${run}"
    
    R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 sim_name='run${run}_' m=5000 disp_prior_df=c(0,c(1,3,5)/10000,c(1,3,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,3,5,10) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/dispersion_error_moderation_real_run.R $ROUT/dispersion_error_moderation_real_run_${data_name}_${count_method}_n${n}.Rout
    
  done

done
done
done



##############################################################################
### Test
##############################################################################



for n in 3
do
  
  echo "n${n}"
  
  R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=10 sim_name='test_' r=5 m=100 disp_prior_df=c(0,c(1,3,5)/10000,c(1,3,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,3,5,10) param_pi_path='$DMPARAMS/prop_kim_kallisto.txt' param_gamma_path='$DMPARAMS/disp_genewise_kim_kallisto_lognormal.txt' param_nm_path='$DMPARAMS/nm_kim_kallisto_lognormal.txt' param_nd_path='$DMPARAMS/nd_common_kim_kallisto.txt'" $RCODE/dispersion_error_moderation_real_run.R $ROUT/dispersion_error_moderation_real_run_n${n}.Rout


done


##############################################################################
### Individual runs
##############################################################################











