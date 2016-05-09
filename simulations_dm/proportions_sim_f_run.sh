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

out_suffix='proportions_decay'
workers=1
disp='disp_common_kim_kallisto'

for n in 3 10
do

  for nm in 1000
  do
    
    
    for run in {1..50}
    do
      
    echo "n${n}_nm${nm}_${run}"

      R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' nr_features=c(3:10) save=FALSE out_suffix='${out_suffix}'" $RCODE/proportions_sim_f_run.R $ROUT/proportions_sim_f_run_n${n}_nm${nm}_${out_suffix}.Rout
      
    done
  done
done



out_suffix='proportions_uniform'
workers=1
disp='disp_common_kim_kallisto'

for n in 3 10
do

  for nm in 1000
  do
    
    
    for run in {1..50}
    do
      
    echo "n${n}_nm${nm}_${run}"

      R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' nr_features=c(3:10) save=FALSE out_suffix='${out_suffix}'" $RCODE/proportions_sim_f_run.R $ROUT/proportions_sim_f_run_n${n}_nm${nm}_${out_suffix}.Rout
      
    done
  done
done




######################
### Test
######################


######################
### Individual run
######################



##############################################################################
### Plot
##############################################################################




n="c(3,10)"
nm="c(1000)"
nd=0
disp='disp_common_kim_kallisto'


out_suffix='proportions_decay'


R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 disp='${disp}' out_suffix='${out_suffix}' pdf_width=10 pdf_height=7 out_dir='proportions_f' strip_text_size=22 text_size=22" $RCODE/proportions_sim_f_plots_run.R $ROUT/proportions_sim_f_plots_run.Rout

tail $ROUT/proportions_sim_f_plots_run.Rout

out_suffix='proportions_uniform'


R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 disp='${disp}' out_suffix='${out_suffix}' pdf_width=10 pdf_height=7 out_dir='proportions_f' strip_text_size=22 text_size=22" $RCODE/proportions_sim_f_plots_run.R $ROUT/proportions_sim_f_plots_run.Rout
   
tail $ROUT/proportions_sim_f_plots_run.Rout



















