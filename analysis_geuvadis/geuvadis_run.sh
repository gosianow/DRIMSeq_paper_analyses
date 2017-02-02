#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_code/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

mkdir $ROUT

###################################################
## Run "prepare_data.R"
###################################################



###################################################
## Run DRIMSeq analysis
###################################################


### Run with permutations: permuted p-values from all the genes

for population in 'CEU' 'YRI'
do

  for chr in {22..1}
  do

    echo "${population}_${chr}"

    R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' chr=${chr}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_${population}_chr${chr}.Rout

  done
done


######################################################################################################
## Run sQTLSeekeR analysis on drimseq counts
######################################################################################################

for population in 'CEU' 'YRI'
do

  R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=1 population='${population}' drimseq_results_path='drimseq_0_3_3_analysis_permutations_all_genes'" $RCODE/geuvadis_sqtlseeker_2_1_drimseq_counts.R $ROUT/geuvadis_sqtlseeker_2_1_drimseq_counts_${population}.Rout

  tail $ROUT/geuvadis_sqtlseeker_2_1_drimseq_counts_${population}.Rout

done


###################################################
## Comparison for FDR = 0.05 and sqtlseeker_2_1_analysis_drimseq_counts
###################################################

###### Comparisons for run with permutations: p-values from all the genes v2 - permutated p-values are saved under results$pvalue_perm


method_out='drimseq_0_3_3_analysis_permutations_all_genes'
comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts'
sqtlseeker_results='sqtlseeker_2_1_analysis_drimseq_counts'
FDR=0.05

### Colors
R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/${comparison_out}'" $RCODE/colors.R $ROUT/colors.Rout


for population in 'CEU' 'YRI'
do

  ### Venn diagrams and upset plots with iCOBRA

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' method_out='${method_out}' comparison_out='${comparison_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_comparison_permutations.R $ROUT/geuvadis_drimseq_0_3_3_comparison_permutations_${population}.Rout
  tail $ROUT/geuvadis_drimseq_0_3_3_comparison_permutations_${population}.Rout


  ### For detected sQTLs: Distance to the closest exon, % within exons, mean gene expression, nr transcripts

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' comparison_out='${comparison_out}' FDR=${FDR} workers=1" $RCODE/geuvadis_drimseq_0_3_3_biol_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_biol_comparison_plots_${population}.Rout
  tail $ROUT/geuvadis_drimseq_0_3_3_biol_comparison_plots_${population}.Rout



  ### CAT (concordance at top) plots for top ranked genes

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' method_out='${method_out}' comparison_out='${comparison_out}' Overlaps_function_path='/home/gosia/R/drimseq_code/help_functions/dm_plotOverlaps.R' CAT_function_path='/home/gosia/R/drimseq_code/help_functions/dm_plotCAT.R' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_comparison_plots_${population}.Rout
  tail $ROUT/geuvadis_drimseq_0_3_3_comparison_plots_${population}.Rout


done



##############################
### Positive controls for FDR = 0.05 and sqtlseeker_2_1_analysis_drimseq_counts
##############################

### Prepare files with validated sQTLs or results from other methods

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/prepare_validated_sqtls.R $ROUT/prepare_validated_sqtls.Rout

tail $ROUT/prepare_validated_sqtls.Rout



### Tables with adj p-values for validated sQTLs; Plots of this tables; Plots of expression of validated sQTLs

method_out='drimseq_0_3_3_analysis_permutations_all_genes'
comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts'
positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_drimseq_counts'
sqtlseeker_results='sqtlseeker_2_1_analysis_drimseq_counts'
FDR=0.05


for population in 'CEU' 'YRI'
do

  # GLIMMPS: PCR validated sQTLs

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt' plot_proportions=TRUE plot_tables=TRUE method_out='${method_out}' comparison_out='${comparison_out}' positive_controls_out='${positive_controls_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr_${population}.Rout

  tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr_${population}.Rout



  # GLIMMPS: GWAS SNPs

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt' plot_proportions=TRUE plot_tables=TRUE method_out='${method_out}' comparison_out='${comparison_out}' positive_controls_out='${positive_controls_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_${population}.Rout

  tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_${population}.Rout



  # GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt' plot_proportions=TRUE plot_tables=TRUE method_out='${method_out}' comparison_out='${comparison_out}' positive_controls_out='${positive_controls_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps_${population}.Rout

  tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps_${population}.Rout



  # GLIMMPS: sQTLs detected by GLIMMPS (~140)

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' valid_path='data/validation/glimmps/glimmps_valid_glimmps.txt' plot_proportions=FALSE plot_tables=FALSE method_out='${method_out}' comparison_out='${comparison_out}' positive_controls_out='${positive_controls_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps_${population}.Rout

  tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps_${population}.Rout



  # GEUVADIS: trQTLs detected within GEUVADIS project for the EUR population all sQTLs

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' valid_path='data/validation/geuvadis/geuvadis_valid_geuvadis_all_EUR.txt' plot_proportions=FALSE plot_tables=FALSE method_out='${method_out}' comparison_out='${comparison_out}' positive_controls_out='${positive_controls_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis_all_EUR_${population}.Rout

  tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis_all_EUR_${population}.Rout


  # GEUVADIS: trQTLs detected within GEUVADIS project for the YRI population all sQTLs

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' valid_path='data/validation/geuvadis/geuvadis_valid_geuvadis_all_YRI.txt' plot_proportions=FALSE plot_tables=FALSE method_out='${method_out}' comparison_out='${comparison_out}' positive_controls_out='${positive_controls_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis_all_YRI_${population}.Rout

  tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis_all_YRI_${population}.Rout


  # GEUVADIS: trQTLs detected within GEUVADIS project for the EUR population best sQTLs

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' valid_path='data/validation/geuvadis/geuvadis_valid_geuvadis_best_EUR.txt' plot_proportions=FALSE plot_tables=FALSE method_out='${method_out}' comparison_out='${comparison_out}' positive_controls_out='${positive_controls_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis_best_EUR_${population}.Rout

  tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis_best_EUR_${population}.Rout


  # GEUVADIS: trQTLs detected within GEUVADIS project for the YRI population best sQTLs

  R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='${population}' valid_path='data/validation/geuvadis/geuvadis_valid_geuvadis_best_YRI.txt' plot_proportions=FALSE plot_tables=FALSE method_out='${method_out}' comparison_out='${comparison_out}' positive_controls_out='${positive_controls_out}' sqtlseeker_results='${sqtlseeker_results}' FDR=${FDR}" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis_best_YRI_${population}.Rout

  tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis_best_YRI_${population}.Rout



done




### Plots the structure of the validated genes with gviz

# GLIMMPS: PCR validated sQTLs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt' positive_controls_out='${positive_controls_out}'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_gviz.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_pcr.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_pcr.Rout

# GLIMMPS: GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt'  positive_controls_out='${positive_controls_out}'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_gviz.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_gwas.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_gwas.Rout


# GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt' positive_controls_out='${positive_controls_out}'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_gviz.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_gwas_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_gwas_glimmps.Rout




############################################################################################################



############################################################################################################
