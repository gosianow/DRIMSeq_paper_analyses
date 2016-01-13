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

### min_feature_expr=c(0,5,10,20,50,100) min_feature_prop=NULL

workers=1

for n in 3
do
  
for data_name in 'kim' 'brooks'
do 
  
for count_method in 'kallisto' 'htseq'
do
    
    for run in {1..50}
    do
      
      echo "${data_name}_${count_method}_n${n}_${run}"
      
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=5000 n=${n} min_feature_expr=c(0,5,10,20,50,100) min_feature_prop=NULL param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}_fcutoff.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/filtering_real_run.R $ROUT/filtering_real_run_${data_name}_${count_method}_n${n}.Rout
    
  done

done
done
done



### min_feature_expr=NULL min_feature_prop=c(0,0.001,0.005,0.01,0.05)

workers=1

for n in 3
do
  
for data_name in 'kim' 'brooks'
do 
  
for count_method in 'kallisto' 'htseq'
do
    
    for run in {101..150}
    do
      
      echo "${data_name}_${count_method}_n${n}_${run}"
      
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=5000 n=${n} min_feature_expr=NULL min_feature_prop=c(0,0.001,0.005,0.01,0.05) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}_fcutoff.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/filtering_real_run.R $ROUT/filtering_real_run_${data_name}_${count_method}_n${n}_prop.Rout
    
  done

done
done
done



################################
### Test
################################

# workers=5

# for n in 3
# do
  
# for data_name in 'brooks'
# do 
  
# for count_method in 'kallisto'
# do
    
#     for run in {1..1}
#     do
      
#       echo "${data_name}_${count_method}_n${n}_${run}"
      
    
#     R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='test_' run='run${run}' m=500 n=${n} min_feature_expr=c(0,20) min_feature_prop=rep(0,2) param_pi_path='$DMPARAMS/${data_name}_${count_method}/prop_${data_name}_${count_method}_fcutoff.txt' param_gamma_path='$DMPARAMS/${data_name}_${count_method}/disp_genewise_${data_name}_${count_method}_lognormal.txt' param_nm_path='$DMPARAMS/${data_name}_${count_method}/nm_${data_name}_${count_method}_lognormal.txt' param_nd_path='$DMPARAMS/${data_name}_${count_method}/nd_common_${data_name}_${count_method}.txt'" $RCODE/filtering_real_run.R $ROUT/filtering_real_run_${data_name}_${count_method}_n${n}.Rout
    
#   done

# done
# done
# done


################################
### Individual runs
################################


##############################################################################
### Plot
##############################################################################



n=3
nm="c('nm_brooks_kallisto_lognormal','nm_brooks_htseq_lognormal','nm_kim_kallisto_lognormal','nm_kim_htseq_lognormal')"
nd="c('nd_common_brooks_kallisto','nd_common_brooks_htseq','nd_common_kim_kallisto','nd_common_kim_htseq')"
prop="c('prop_brooks_kallisto_fcutoff','prop_brooks_htseq_fcutoff','prop_kim_kallisto_fcutoff','prop_kim_htseq_fcutoff')"
disp="c('disp_genewise_brooks_kallisto_lognormal','disp_genewise_brooks_htseq_lognormal','disp_genewise_kim_kallisto_lognormal','disp_genewise_kim_htseq_lognormal')"
data_name="c('brooks','kim')"
count_method="c('kallisto','htseq')"


out_suffix='filtering_real_min_feature_prop'


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=${nd} prop=${prop} disp=${disp} data_name=${data_name} count_method=${count_method} out_suffix='${out_suffix}'" $RCODE/filtering_real_plots_run.R $ROUT/filtering_real_plots_run.Rout
     
    



out_suffix='filtering_real_min_feature_expr'


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=${nd} prop=${prop} disp=${disp} data_name=${data_name} count_method=${count_method} out_suffix='${out_suffix}'" $RCODE/filtering_real_plots_run.R $ROUT/filtering_real_plots_run.Rout
     
  









