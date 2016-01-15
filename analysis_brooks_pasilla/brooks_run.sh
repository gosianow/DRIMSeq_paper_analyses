#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_brooks_pasilla
RWD=/home/Shared/data/seq/brooks_pasilla
ROUT=$RWD/Rout
ANNOTATION=/home/Shared/data/annotation/Drosophila/Ensembl70

mkdir $ROUT

## Run R scripts


###############################################################################
### Download data + Mapping with tophat2
###############################################################################

### ! First save SraRunInfo.csv file in "$RWD/3_metadata" directory !


R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' DNA_fasta='$ANNOTATION/genome/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa'" $RCODE/brooks_download.R $ROUT/brooks_download.Rout



###############################################################################
### Kallisto
###############################################################################


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf_path='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' cDNA_fasta='$ANNOTATION/cDNA/Drosophila_melanogaster.BDGP5.70.cdna.all.fa'" $RCODE/brooks_kallisto.R $ROUT/brooks_kallisto.Rout



R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf_path='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf'" $RCODE/brooks_kallisto_filter_gtf.R $ROUT/brooks_kallisto_filter_gtf.Rout


###############################################################################
### DEXSeq
###############################################################################


### Index BAM files

for i in 'GSM461176' 'GSM461177' 'GSM461178' 'GSM461179' 'GSM461180' 'GSM461181' 'GSM461182'
 do 
  samtools index $RWD/1_reads/tophat_2.0.14/${i}/accepted_hits.bam
done


### Run HTSeq

R214 CMD BATCH --no-save --no-restore "--args gtf='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' count_method='htseq'" $RCODE/brooks_htseq.R $ROUT/brooks_htseq_htseq.Rout


R214 CMD BATCH --no-save --no-restore "--args gtf='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70_kallistoest_atleast5.gtf' count_method='htseqprefiltered5'" $RCODE/brooks_htseq.R $ROUT/brooks_htseq_htseqprefiltered5.Rout


### Run DEXSeq

for model in 'model_full' 'model_full_glm' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R214 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=10 count_method='${count_method}' model='${model}'" $RCODE/brooks_dexseq_run.R $ROUT/brooks_dexseq_run_${model}_${count_method}.Rout

  done
done



###############################################################################
### DRIMSeq
###############################################################################


for model in 'model_full' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do
    
    echo "${model}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=TRUE results_common=FALSE disp_mode_list='grid' disp_moderation_list='none'" $RCODE/brooks_drimseq_0_3_3_run.R $ROUT/brooks_drimseq_0_3_3_run_${model}_${count_method}_grid_none.Rout
    
    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 count_method='${count_method}' model='${model}' dispersion_common=FALSE results_common=FALSE disp_mode_list='grid' disp_moderation_list='common'" $RCODE/brooks_drimseq_0_3_3_run.R $ROUT/brooks_drimseq_0_3_3_run_${model}_${count_method}_grid_common.Rout

  done
done

##############################
### Colors
##############################

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/drimseq_0_3_3_comparison'" $RCODE/colors.R $ROUT/colors.Rout


##############################
### DRIMSeq comparison
##############################


for model in 'model_full' 'model_full_glm' 'model_full_paired' 'model_null1' 'model_null2' 'model_null3'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${model}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}' model='${model}'" $RCODE/brooks_drimseq_0_3_3_comparison_run.R $ROUT/brooks_drimseq_0_3_3_comparison_run.Rout

  done
done


R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/brooks_drimseq_0_3_3_summary.R $ROUT/brooks_drimseq_0_3_3_summary.Rout


##############################
### DRIMSeq positive controls
##############################

### ! Save brooks_validated_genes.txt in "$RWD/5_validation" directory !



for model in 'model_full' 'model_full_glm' 'model_full_paired'
do 
  for count_method in 'kallisto' 'htseq' 'kallistofiltered5' 'htseqprefiltered5'
  do 
  
    echo "${model}_${count_method}"

    R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' gtf_path='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' count_method='${count_method}' model='${model}'" $RCODE/brooks_drimseq_0_3_3_positive_controls.R $ROUT/brooks_drimseq_0_3_3_positive_controls_run.Rout

  done
done





R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' path_gtf='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70.gtf' path_gtf_filtered='$ANNOTATION/gtf/Drosophila_melanogaster.BDGP5.70_kallistoest_atleast5.gtf'" $RCODE/brooks_drimseq_0_3_3_positive_controls_summary.R $ROUT/brooks_drimseq_0_3_3_positive_controls_summary.Rout














