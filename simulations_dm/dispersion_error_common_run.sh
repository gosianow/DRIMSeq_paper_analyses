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

  for nm in 100 1000
  do
    
    for prop in 'prop_q3_uniform' 'prop_q10_uniform' 'prop_q3_kim_kallisto_overall' 'prop_q10_kim_kallisto_overall'
    do 
    
    echo "n${n}_nm${nm}_${prop}"

      R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 sim_name='' r=50 m=100 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/${prop}.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_common_run.R $ROUT/dispersion_error_common_run_n${n}_nm${nm}_${prop}.Rout

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
    
    for prop in 'prop_q3_uniform'
    do 
    
    echo "n${n}_nm${nm}_${prop}"

      R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 sim_name='test_' r=2 m=100 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/${prop}.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_common_run.R $ROUT/dispersion_error_common_run_n${n}_nm${nm}_${prop}.Rout

    done
  done
done


#############################################################################
### Individual run
#############################################################################


for n in 6
do

  for nm in 100 1000
  do
    
    for prop in 'prop_q3_uniform' 'prop_q10_uniform' 'prop_q3_kim_kallisto_overall' 'prop_q10_kim_kallisto_overall'
    do 
    
    echo "n${n}_nm${nm}_${prop}"

      R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 sim_name='' r=50 m=100 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/${prop}.txt' param_gamma_path='$DMPARAMS/disp_common_kim_kallisto.txt'" $RCODE/dispersion_error_common_run.R $ROUT/dispersion_error_common_run_n${n}_nm${nm}_${prop}.Rout

    done
  done
done




