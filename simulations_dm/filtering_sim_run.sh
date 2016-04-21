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
disp='disp_common_kim_kallisto'
prop='prop_q20_kim_kallisto_fcutoff'

for n in 3 6
do

  for nm in 10000 100000
  do
    
    
    for run in {1..50}
    do
      
    echo "n${n}_nm${nm}_${prop}_${run}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' max_features=c(Inf,18,16,14,12,10,8,5)" $RCODE/filtering_sim_run.R $ROUT/filtering_sim_run_n${n}_nm${nm}_${prop}.Rout
      
    done
  done
done




######################
### Test
######################

# workers=5
# disp='disp_common_kim_kallisto'
# prop='prop_q20_kim_kallisto_fcutoff'

# for n in 3
# do

#   for nm in 1000
#   do
    
    
#     for run in {1..1}
#     do
      
#     echo "n${n}_nm${nm}_${prop}_${run}"

#       R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='test_' run='run${run}' m=100 n=${n} nm=${nm} nd=0 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' max_features=c(Inf,18)" $RCODE/filtering_sim_run.R $ROUT/filtering_sim_run_n${n}_nm${nm}_${prop}.Rout
      
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
nm="c(10000,100000)"
nd=0
prop='prop_q20_kim_kallisto_fcutoff'
disp='disp_common_kim_kallisto'
out_suffix='filtering'


R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop='${prop}' disp='${disp}' out_suffix='${out_suffix}' pdf_width=7 pdf_height=7" $RCODE/filtering_sim_plots_run.R $ROUT/filtering_sim_plots_run.Rout
     
tail $ROUT/filtering_sim_plots_run.Rout







