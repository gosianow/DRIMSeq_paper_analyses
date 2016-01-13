#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters_drimseq_0_3_3

mkdir $ROUT


##############################################################################
### Run
##############################################################################

workers=1

for n in 3 6
do
  
for data_name in 'kim' 'brooks'
do 
  
for count_method in 'kallisto' 'htseq'
do
    
    for run in {1..25}
    do
      
      echo "${data_name}_${count_method}_n${n}_${run}"
      
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=5000 n=${n} disp_prior_df=c(0,c(1,5)/10000,c(1,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,5,10) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}_fcutoff.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/moderation_real_run.R $ROUT/moderation_real_run_${data_name}_${count_method}_n${n}.Rout
    
  done

done
done
done



################################
### Test
################################

# RCODE=/home/gosia/R/drimseq_paper/simulations_dm
# RWD=/home/gosia/multinomial_project/simulations_dm/drimseq
# ROUT=$RWD/Rout
# DMPARAMS=$RWD/dm_parameters_drimseq_0_3_3


# workers=10

# for n in 3
# do
  
# for data_name in 'brooks'
# do 
  
# for count_method in 'kallisto'
# do
    
#     for run in {1..1}
#     do
      
#       echo "${data_name}_${count_method}_n${n}_${run}"
      
    
#     R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='test_' run='run${run}' m=5000 n=${n} disp_prior_df=c(0,c(1,5)/10000,c(1,5)/1000,c(1,3,5)/100,c(1,3,5)/10,1,5,10) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}_fcutoff.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/moderation_real_run.R $ROUT/moderation_real_run_${data_name}_${count_method}_n${n}.Rout
    
#   done

# done
# done
# done


################################
### Individual runs
################################


##############################################################################
### Plots
##############################################################################

n="c(3,6)"
nm="c('nm_brooks_kallisto_lognormal','nm_brooks_htseq_lognormal','nm_kim_kallisto_lognormal','nm_kim_htseq_lognormal')"
nd="c('nd_common_brooks_kallisto','nd_common_brooks_htseq','nd_common_kim_kallisto','nd_common_kim_htseq')"
prop="c('prop_brooks_kallisto_fcutoff','prop_brooks_htseq_fcutoff','prop_kim_kallisto_fcutoff','prop_kim_htseq_fcutoff')"
disp="c('disp_genewise_brooks_kallisto_lognormal','disp_genewise_brooks_htseq_lognormal','disp_genewise_kim_kallisto_lognormal','disp_genewise_kim_htseq_lognormal')"
data_name="c('brooks','kim')"
count_method="c('kallisto','htseq')"
out_suffix='moderation_real'


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=${nd} prop=${prop} disp=${disp} data_name=${data_name} count_method=${count_method} out_suffix='${out_suffix}' pdf_width=14 pdf_height=7" $RCODE/moderation_real_plots_run.R $ROUT/moderation_real_plots_run_n${n}.Rout
    

















