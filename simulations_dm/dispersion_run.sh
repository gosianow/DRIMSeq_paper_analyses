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

disp='disp_genewise_kim_kallisto_lognormal'
workers=1


for n in 3 6
do

  for nm in 100 1000
  do
    
    for prop in 'prop_q3_uniform' 'prop_q10_uniform' 'prop_q3_kim_kallisto_overall' 'prop_q10_kim_kallisto_overall'
    do 
    
    for run in {1..25}
    do
      
    echo "n${n}_nm${nm}_${prop}_${run}"
    
      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt'" $RCODE/dispersion_run.R $ROUT/dispersion_run_${disp}_n${n}_nm${nm}_${prop}.Rout
      
    done

    done
  done
done



disp='disp_common_kim_kallisto'
workers=1


for n in 3 6
do

  for nm in 100 1000
  do
    
    for prop in 'prop_q3_uniform' 'prop_q10_uniform' 'prop_q3_kim_kallisto_overall' 'prop_q10_kim_kallisto_overall'
    do 
    
    for run in {1..25}
    do
      
    echo "n${n}_nm${nm}_${prop}_${run}"
    
      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt'" $RCODE/dispersion_run.R $ROUT/dispersion_run_${disp}_n${n}_nm${nm}_${prop}.Rout
      
    done

    done
  done
done






################################
### Test
################################

disp='disp_common_kim_kallisto'
workers=4


for n in 3
do

  for nm in 100
  do
    
    for prop in 'prop_q3_uniform'
    do 
    
    for run in {1..3}
    do
      
    echo "n${n}_nm${nm}_${prop}_${run}"
    
      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='test_' run='run${run}' m=100 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt'" $RCODE/dispersion_run.R $ROUT/dispersion_run_${disp}_n${n}_nm${nm}_${prop}.Rout
      
    done

    done
  done
done





##############################################################################
### Plot
##############################################################################


n="c(3,6)"
nm="c(100,1000)"
nd=0
prop="c('prop_q3_uniform','prop_q10_uniform','prop_q3_kim_kallisto_overall','prop_q10_kim_kallisto_overall')"
disp="c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')"


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop=${prop} disp=${disp}" $RCODE/dispersion_plots_run.R $ROUT/dispersion_plots_run.Rout
































