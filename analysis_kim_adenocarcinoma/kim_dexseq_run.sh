## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout

# mkdir $ROUT

## Run R scripts


for i in 'kallisto' 'htseq'
do 
  for j in 'model_null_normal2' 'model_null_tumor2'
  do

R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${i}' model='${j}'" $RCODE/kim_dexseq_run.R $ROUT/kim_dexseq_run_${i}_${j}.Rout

  done
done







