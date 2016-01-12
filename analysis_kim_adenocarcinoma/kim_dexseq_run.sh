## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_kim_adenocarcinoma
RWD=/home/Shared/data/seq/kim_adenocarcinoma
ROUT=$RWD/Rout

# mkdir $ROUT

## Run R scripts


### Index BAM files

for i in 'GSM927308' 'GSM927310' 'GSM927312' 'GSM927314' 'GSM927316' 'GSM927318' 'GSM927309' 'GSM927311' 'GSM927313' 'GSM927315' 'GSM927317' 'GSM927319'
 do 
  samtools index $RWD/1_reads/tophat_insilicodb/${i}/accepted_hits.bam
done


### Run HTSeq


R214 CMD BATCH --no-save --no-restore "--args gtf='/home/Shared_penticton/data/annotation/Human/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71.gtf' count_method='htseq'" $RCODE/kim_htseq.R $ROUT/kim_htseq_htseq.Rout


R214 CMD BATCH --no-save --no-restore "--args gtf='/home/Shared_penticton/data/annotation/Human/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71_kallistoest_atleast5.gtf' count_method='htseqprefiltered5'" $RCODE/kim_htseq.R $ROUT/kim_htseq_htseqprefiltered5.Rout




### Run DEXSeq

for model in 'model_full' 'model_full_glm' 'model_null_normal1' 'model_null_normal2' 'model_null_tumor1' 'model_null_tumor2'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do

    echo "${model}_${count_method}"

    R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}'" $RCODE/kim_dexseq_run.R $ROUT/kim_dexseq_run_${model}_${count_method}.Rout

  done
done


###############################################################################
### Individual runs
###############################################################################



for model in 'model_full_glm'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do

    echo "${model}_${count_method}"

    R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}'" $RCODE/kim_dexseq_run.R $ROUT/kim_dexseq_run_${model}_${count_method}.Rout

  done
done
















