#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_code/simulations_dm
RWD=/home/gosia/multinomial_project/simulations_dm/drimseq
ROUT=$RWD/Rout
DMPARAMS=$RWD/dm_parameters_drimseq_0_3_3

mkdir $ROUT

##############################################################################
### Run DRIMSeq version with fixed priorN
##############################################################################

disp='disp_genewise_kim_kallisto_lognormal'
workers=1


for n in 3 6 12
do

  for nm in 100 1000
  do

    for prop in 'prop_q3_uniform' 'prop_q10_uniform' 'prop_q3_kim_kallisto_fcutoff' 'prop_q10_kim_kallisto_fcutoff'
    do

    for run in {1..50}
    do

    echo "n${n}_nm${nm}_${prop}_${run}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' out_dir='core_fp' out_suffix='core_fp'" $RCODE/core_fp_run.R $ROUT/core_fp_run_${disp}_n${n}_nm${nm}_${prop}.Rout

    done

    done
  done
done



disp='disp_common_kim_kallisto'
workers=1


for n in 3 6 12
do

  for nm in 100 1000
  do

    for prop in 'prop_q3_uniform' 'prop_q10_uniform' 'prop_q3_kim_kallisto_fcutoff' 'prop_q10_kim_kallisto_fcutoff'
    do

    for run in {1..50}
    do

    echo "n${n}_nm${nm}_${prop}_${run}"

      R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' out_dir='core_fp' out_suffix='core_fp'" $RCODE/core_fp_run.R $ROUT/core_fp_run_${disp}_n${n}_nm${nm}_${prop}.Rout

    done

    done
  done
done




##############################################################################
### Run DRIMSeq version with auto calculated priorN - equal median for common moderation
##############################################################################

disp='disp_genewise_kim_kallisto_lognormal'
workers=1


for n in 3 6 12
do

  for nm in 100 1000
  do

    for prop in 'prop_q3_uniform' 'prop_q10_uniform' 'prop_q3_kim_kallisto_fcutoff' 'prop_q10_kim_kallisto_fcutoff'
    do

    for run in {1..50}
    do

    echo "n${n}_nm${nm}_${prop}_${run}"

      R32devloc CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' out_dir='core_fp_auto_moderation' out_suffix='core_fp'" $RCODE/core_fp_run.R $ROUT/core_fp_run_${disp}_n${n}_nm${nm}_${prop}.Rout

    done

    done
  done
done



disp='disp_common_kim_kallisto'
workers=1


for n in 3 6 12
do

  for nm in 100 1000
  do

    for prop in 'prop_q3_uniform' 'prop_q10_uniform' 'prop_q3_kim_kallisto_fcutoff' 'prop_q10_kim_kallisto_fcutoff'
    do

    for run in {1..50}
    do

    echo "n${n}_nm${nm}_${prop}_${run}"

      R32devloc CMD BATCH --no-save --no-restore "--args rwd='$RWD' simulation_script='$RCODE/dm_simulate.R' workers=${workers} sim_name='' run='run${run}' m=1000 n=${n} nm=${nm} nd=0 disp_prior_df=0.1 param_pi_path='$DMPARAMS/kim_kallisto/${prop}.txt' param_gamma_path='$DMPARAMS/kim_kallisto/${disp}.txt' out_dir='core_fp_auto_moderation' out_suffix='core_fp'" $RCODE/core_fp_run.R $ROUT/core_fp_run_${disp}_n${n}_nm${nm}_${prop}.Rout

    done

    done
  done
done



################################
### Test
################################

################################
### Individual runs
################################

##############################################################################
### Plots for core_fp
##############################################################################

n="c(3,6,12)"
nm="c(100,1000)"
nd=0
prop="c('prop_q3_uniform','prop_q3_kim_kallisto_fcutoff','prop_q10_uniform','prop_q10_kim_kallisto_fcutoff')"
disp="c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')"


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop=${prop} disp=${disp} pdf_width=20 pdf_height=10 fig_name='all_' out_dir='core_fp' out_suffix='core_fp'" $RCODE/core_fp_plots_run.R $ROUT/core_fp_plots_run.Rout


n="c(3,12)"
nm="c(100,1000)"
nd=0
prop="c('prop_q3_uniform','prop_q3_kim_kallisto_fcutoff','prop_q10_uniform','prop_q10_kim_kallisto_fcutoff')"
disp="c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')"


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop=${prop} disp=${disp} pdf_width=14 pdf_height=8 fig_name='3_12_100_1000_' out_dir='core_fp' out_suffix='core_fp'" $RCODE/core_fp_plots_run.R $ROUT/core_fp_plots_run.Rout


n="c(3,12)"
nm="c(1000)"
nd=0
prop="c('prop_q3_uniform','prop_q3_kim_kallisto_fcutoff','prop_q10_uniform','prop_q10_kim_kallisto_fcutoff')"
disp="c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')"


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop=${prop} disp=${disp} pdf_width=8 pdf_height=10 fig_name='3_12_1000_' out_dir='core_fp' out_suffix='core_fp'" $RCODE/core_fp_plots_run.R $ROUT/core_fp_plots_run.Rout



n="c(6)"
nm="c(1000)"
nd=0
prop="c('prop_q3_kim_kallisto_fcutoff','prop_q10_kim_kallisto_fcutoff')"
disp="c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')"


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop=${prop} disp=${disp} pdf_width=7 pdf_height=7 fig_name='6_1000_' out_dir='core_fp' out_suffix='core_fp'" $RCODE/core_fp_plots_run.R $ROUT/core_fp_plots_run.Rout



n="c(3)"
nm="c(1000)"
nd=0
prop="c('prop_q3_kim_kallisto_fcutoff','prop_q10_kim_kallisto_fcutoff')"
disp="c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')"


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop=${prop} disp=${disp} pdf_width=7 pdf_height=7 fig_name='3_1000_' out_dir='core_fp' out_suffix='core_fp'" $RCODE/core_fp_plots_run.R $ROUT/core_fp_plots_run.Rout





##############################################################################
### Plots for core_fp_auto_moderation
##############################################################################

n="c(3,12)"
nm="c(1000)"
nd=0
prop="c('prop_q3_uniform','prop_q3_kim_kallisto_fcutoff','prop_q10_uniform','prop_q10_kim_kallisto_fcutoff')"
disp="c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')"


R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop=${prop} disp=${disp} pdf_width=8 pdf_height=10 fig_name='3_12_1000_' out_dir='core_fp_auto_moderation' out_suffix='core_fp' strip_text_size=10 text_size=16" $RCODE/core_fp_plots_run.R $ROUT/core_fp_plots_run.Rout



n="c(3)"
nm="c(1000)"
nd=0
prop="c('prop_q3_kim_kallisto_fcutoff','prop_q10_kim_kallisto_fcutoff')"
disp="c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')"


R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' sim_name='' n=${n} nm=${nm} nd=0 prop=${prop} disp=${disp} pdf_width=7 pdf_height=7 fig_name='3_1000_' out_dir='core_fp_auto_moderation' out_suffix='core_fp' strip_text_size=14 text_size=16" $RCODE/core_fp_plots_run.R $ROUT/core_fp_plots_run.Rout

tail $ROUT/core_fp_plots_run.Rout














#
