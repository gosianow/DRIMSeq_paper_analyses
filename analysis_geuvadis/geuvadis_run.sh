#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

mkdir $ROUT

###################################################
## Run "prepare_data.R"
###################################################


###################################################
## Run DRIMSeq analysis
###################################################


for n in {22..1}
do 

echo "${n}"

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 population='CEU' chr=${n}" $RCODE/geuvadis_drimseq_0_3_3_run.R $ROUT/geuvadis_drimseq_0_3_3_run_CEU_chr${n}.Rout

done


###################################################
## Run DRIMSeq analysis using F test - R32dev
###################################################


for n in 22
do 

echo "${n}"

R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=6 population='CEU' chr=${n}" $RCODE/geuvadis_drimseq_0_3_3_run_f.R $ROUT/geuvadis_drimseq_0_3_3_run_f_CEU_chr${n}.Rout

done


###################################################
## Run DRIMSeq analysis on random SNP-gene pairs 
###################################################


for n in {22..1}
do 

echo "${n}"

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 population='CEU' chr=${n}" $RCODE/geuvadis_drimseq_0_3_3_random.R $ROUT/geuvadis_drimseq_0_3_3_random_CEU_chr${n}.Rout

done


###################################################
## Run sQTLSeekeR analysis
###################################################


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4" $RCODE/geuvadis_sqtlseeker_2_1_run.R $ROUT/geuvadis_sqtlseeker_2_1_run.Rout



###################################################
## Run sQTLSeekeR analysis on random SNP-gene pairs
###################################################



R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4" $RCODE/geuvadis_sqtlseeker_2_1_random.R $ROUT/geuvadis_sqtlseeker_2_1_random.Rout
























