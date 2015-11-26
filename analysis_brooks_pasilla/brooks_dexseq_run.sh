## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_brooks_pasilla
RWD=/home/Shared/data/seq/brooks_pasilla
ROUT=$RWD/Rout

mkdir $ROUT

## Run R scripts


### Index BAM files

for i in 'GSM461176' 'GSM461177' 'GSM461178' 'GSM461179' 'GSM461180' 'GSM461181' 'GSM461182'
 do 
  samtools index $RWD/1_reads/tophat_2.0.9/${i}/accepted_hits.bam
done


### Run HTSeq

R214 CMD BATCH --no-save --no-restore "--args gtf='/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70_kallistoest_atleast5.gtf' count_method='htseqprefiltered5'" $RCODE/brooks_htseq.R $ROUT/brooks_htseq_htseqprefiltered5.Rout


### Run DEXSeq

for model in 'model_full' 'model_full_glm' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}'" $RCODE/brooks_dexseq_run.R $ROUT/brooks_dexseq_run_${model}_${count_method}.Rout

  done
done



###############################################################################
### Individual runs
###############################################################################


for model in 'model_full' 'model_full_glm' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
do 
  for count_method in 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=2 count_method='${count_method}' model='${model}'" $RCODE/brooks_dexseq_run.R $ROUT/brooks_dexseq_run_${model}_${count_method}.Rout

  done
done
























