## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_brooks_pasilla
RWD=/home/Shared/data/seq/brooks_pasilla
ROUT=$RWD/Rout

mkdir $ROUT

## Run R scripts

for i in 'kallisto' 'htseq'
do 
  for j in 'model_full' 'model_full_glm' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
  do

R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${i}' model='${j}'" $RCODE/brooks_dexseq_run.R $ROUT/brooks_dexseq_run_${i}_${j}.Rout

  done
done



