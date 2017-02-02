#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_code/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters_drimseq_0_3_3

mkdir $ROUT

##############################################################################
### Run
##############################################################################

out_suffix='proportions_decay'
workers=1
disp='disp_genewise_kim_kallisto_lognormal'

for n in 3 6
do

  for nm in 1000 100000
  do


    for run in {1..50}
    do

    echo "n${n}_nm${nm}_${run}_${out_suffix}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' nr_features=c(2,4,6,8,10,13,16,19) out_suffix='${out_suffix}'" $RCODE/proportions_and_moderation_sim_run.R $ROUT/proportions_and_moderation_sim_run_n${n}_nm${nm}_${out_suffix}.Rout

    done
  done
done




out_suffix='proportions_uniform'
workers=1
disp='disp_genewise_kim_kallisto_lognormal'

for n in 3 6
do

  for nm in 1000 100000
  do


    for run in {1..50}
    do

    echo "n${n}_nm${nm}_${run}_${out_suffix}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' nr_features=c(2,4,6,8,10,13,16,19) out_suffix='${out_suffix}'" $RCODE/proportions_and_moderation_sim_run.R $ROUT/proportions_and_moderation_sim_run_n${n}_nm${nm}_${out_suffix}.Rout

    done
  done
done


######################
### Test
######################

# out_suffix='proportions_uniform'
# workers=5
# disp='disp_genewise_kim_kallisto_lognormal'

# for n in 3
# do

#   for nm in 1000
#   do


#     for run in {1..2}
#     do

#     echo "n${n}_nm${nm}_${run}_${out_suffix}"

#       R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='test_' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' nr_features=c(2,10) out_suffix='${out_suffix}'" $RCODE/proportions_and_moderation_sim_run.R $ROUT/proportions_and_moderation_sim_run_n${n}_nm${nm}_${out_suffix}.Rout

#     done
#   done
# done



######################
### Individual run
######################


##############################################################################
### Plot
##############################################################################



n="c(3,6)"
nm="c(1000,100000)"
nd=0
disp='disp_genewise_kim_kallisto_lognormal'


out_suffix='proportions_decay'


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 disp='${disp}' out_suffix='${out_suffix}' pdf_width=7 pdf_height=7" $RCODE/proportions_and_moderation_sim_plots_run.R $ROUT/proportions_and_moderation_sim_plots_run.Rout


out_suffix='proportions_uniform'


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 disp='${disp}' out_suffix='${out_suffix}' pdf_width=7 pdf_height=7" $RCODE/proportions_and_moderation_sim_plots_run.R $ROUT/proportions_and_moderation_sim_plots_run.Rout
