#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters

# mkdir $ROUT

## Run R scripts

# common dispersion from kallisto


for n in 3 6
do

  for nm in 1000
  do
    
    for prop in 'prop_q15_kim_kallisto_overall'
    do 
    
    for run in {1..50}
    do
      
    echo "n${n}_nm${nm}_${prop}_${run}"

      R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='run${run}_' m=1000 n=${n} nm=${nm} nd=0 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/disp_common_kim_kallisto.txt' max_features=c(Inf,13,12,10,8)" $RCODE/dispersion_error_filtering_run.R $ROUT/dispersion_error_filtering_run_n${n}_nm${nm}_${prop}.Rout
      
    done
    done
  done
done




for n in 3
do

  for nm in 1000
  do
    
    for prop in 'prop_q20_kim_kallisto_fcutoff'
    do 
    
    for run in {1..50}
    do
      
    echo "n${n}_nm${nm}_${prop}_${run}"

      R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 sim_name='run${run}_' m=1000 n=${n} nm=${nm} nd=0 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/disp_common_kim_kallisto.txt' max_features=c(Inf,18,16,14,12,10,8)" $RCODE/dispersion_error_filtering_run.R $ROUT/dispersion_error_filtering_run_n${n}_nm${nm}_${prop}.Rout
      
    done
    done
  done
done



#############################################################################
### Test
#############################################################################

for n in 3
do

  for nm in 100
  do
    
    for prop in 'prop_q15_kim_kallisto_overall'
    do 
    
    for run in {1..3}
    do
      
    echo "n${n}_nm${nm}_${prop}_${run}"

      R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='run${run}_' m=100 n=${n} nm=${nm} nd=0 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/disp_common_kim_kallisto.txt' max_features=c(Inf,13,12,10,8)" $RCODE/dispersion_error_filtering_run.R $ROUT/dispersion_error_filtering_run_n${n}_nm${nm}_${prop}.Rout
      
    done
    done
  done
done


#############################################################################
### Individual run
#############################################################################



RCODE=/home/gosia/R/drimseq_paper/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters


for n in 3
do

  for nm in 1000
  do
    
    for prop in 'prop_q15_kim_kallisto_overall'
    do 
    
    for run in {1..50}
    do
      
    echo "n${n}_nm${nm}_${prop}_${run}"

      R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='run${run}_' m=1000 n=${n} nm=${nm} nd=0 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/disp_common_kim_kallisto.txt' max_features=c(Inf,13,12,10,8)" $RCODE/dispersion_error_filtering_run.R $ROUT/dispersion_error_filtering_run_n${n}_nm${nm}_${prop}.Rout
      
    done
    done
  done
done






















