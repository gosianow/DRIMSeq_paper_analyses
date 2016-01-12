#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters_drimseq_0_3_3

# mkdir $ROUT

##############################################################################
### Run
##############################################################################

workers=1
disp='disp_common_kim_kallisto'

for n in 3 6
do

  for nm in 1000 10000
  do
    
    
    for run in {1..50}
    do
      
    echo "n${n}_nm${nm}_${run}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' nr_features=c(3,5,7,10,13,15,17)" $RCODE/proportions_sim_run.R $ROUT/proportions_sim_run_n${n}_nm${nm}.Rout
      
    done
  done
done




######################
### Test
######################

workers=5
disp='disp_common_kim_kallisto'

for n in 3
do

  for nm in 1000
  do
    
    
    for run in {1..1}
    do
      
    echo "n${n}_nm${nm}_${run}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='test_' run='run${run}' m=100 n=${n} nm=${nm} nd=0 param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' nr_features=c(3,5)" $RCODE/proportions_sim_run.R $ROUT/proportions_sim_run_n${n}_nm${nm}.Rout
      
    done
  done
done



######################
### Individual run
######################

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters_drimseq_0_3_3

workers=1
disp='disp_common_kim_kallisto'

for n in 6
do

  for nm in 1000
  do
    
    
    for run in {1..50}
    do
      
    echo "n${n}_nm${nm}_${run}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' nr_features=c(3,5,7,10,13,15,17)" $RCODE/proportions_sim_run.R $ROUT/proportions_sim_run_n${n}_nm${nm}.Rout
      
    done
  done
done


##############################################################################
### Plot
##############################################################################



n="c(3,6)"
nm="c(10000,1000)"
nd=0
disp='disp_common_kim_kallisto'
out_suffix='proportions_decay'


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 disp='${disp}' out_suffix='${out_suffix}' pdf_width=7 pdf_height=7" $RCODE/proportions_sim_plots_run.R $ROUT/proportions_sim_plots_run.Rout
     

out_suffix='proportions_uniform'


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 disp='${disp}' out_suffix='${out_suffix}' pdf_width=7 pdf_height=7" $RCODE/proportions_sim_plots_run.R $ROUT/proportions_sim_plots_run.Rout
   









