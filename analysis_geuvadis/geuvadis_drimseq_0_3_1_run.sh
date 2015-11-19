## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

mkdir $ROUT

## Run R scripts

for n in {1..22}
do 
R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 population='CEU' chr=${n}" $RCODE/geuvadis_drimseq_0_3_1_run.R $ROUT/geuvadis_drimseq_0_3_1_run_CEU_chr${n}.Rout
done


for n in {1..22}
do 
R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5 population='YRI' chr=${n}" $RCODE/geuvadis_drimseq_0_3_1_run.R $ROUT/geuvadis_drimseq_0_3_1_run_YRI_chr${n}.Rout
done

