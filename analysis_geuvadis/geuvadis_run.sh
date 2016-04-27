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


### Run with permutations: p-values from all the genes

for n in {22..1}
do 

echo "${n}"

R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=10 population='CEU' chr=${n}" $RCODE/geuvadis_drimseq_0_3_3_run_permutations.R $ROUT/geuvadis_drimseq_0_3_3_run_permutations_CEU_chr${n}.Rout

done



###################################################
## Run sQTLSeekeR analysis
###################################################


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4" $RCODE/geuvadis_sqtlseeker_2_1_run.R $ROUT/geuvadis_sqtlseeker_2_1_run.Rout


###################################################
## Comparison for FDR = 0.05
###################################################

###### Comparisons for run with permutations: p-values from all the genes v2 - permutated p-values are saved under results$pvalue_perm



### Colors 
R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/drimseq_0_3_3_comparison_permutations_all_genes'" $RCODE/colors.R $ROUT/colors.Rout


### Venn diagrams and upset plots with iCOBRA

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_comparison_permutations.R $ROUT/geuvadis_drimseq_0_3_3_comparison_permutations.Rout


### For detected sQTLs: Distance to the closest exon, % within exons, mean gene expression, nr transcripts

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes' FDR=0.05 workers=10" $RCODE/geuvadis_drimseq_0_3_3_biol_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_biol_comparison_plots.Rout


### Plots of the overlap versus number of top ranked genes
### CAT (concordance at top) plots

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes' Overlaps_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotOverlaps.R' CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_comparison_plots.Rout


##############################
### Positive controls for FDR = 0.05
##############################

### Prepare lists of PCR validated sQTLs and GWAS SNPs from GLiMMPS paper
### Run the "prepare_validated_sqtls.R" script




### Tables with adj p-values for validated sQTLs; Plots of this tables; Plots of expression of validated sQTLs

# GLIMMPS: PCR validated sQTLs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr.Rout



# GLIMMPS: GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas.Rout



# GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps.Rout



# GLIMMPS: sQTLs detected by GLIMMPS (~140)

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_glimmps.txt' plot_proportions=FALSE plot_tables=FALSE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps.Rout



# GEUVADIS: trQTLs detected within GEUVADIS project for the EUR population

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/geuvadis/geuvadis_valid_geuvadis.txt' plot_proportions=FALSE plot_tables=FALSE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis.Rout






### Plots the structure of the validated genes with gviz

# GLIMMPS: PCR validated sQTLs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_gviz.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_pcr.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_pcr.Rout

# GLIMMPS: GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt'  positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_gviz.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_gwas.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_gwas.Rout


# GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_gviz.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_gwas_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gviz_gwas_glimmps.Rout




###################################################
## Comparison for FDR = 0.1
###################################################

###### Comparisons for run with permutations: p-values from all the genes v2 - permutated p-values are saved under results$pvalue_perm



### Colors 
R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/drimseq_0_3_3_comparison_permutations_all_genes_fdr010'" $RCODE/colors.R $ROUT/colors.Rout


### Venn diagrams and upset plots with iCOBRA

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_fdr010' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.1" $RCODE/geuvadis_drimseq_0_3_3_comparison_permutations.R $ROUT/geuvadis_drimseq_0_3_3_comparison_permutations.Rout



### For detected sQTLs: Distance to the closest exon, % within exons, mean gene expression, nr transcripts

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_fdr010' FDR=0.1 workers=10" $RCODE/geuvadis_drimseq_0_3_3_biol_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_biol_comparison_plots.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_biol_comparison_plots.Rout



### Plots of the overlap versus number of top ranked genes
### CAT (concordance at top) plots

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_fdr010' Overlaps_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotOverlaps.R' CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R' FDR=0.1" $RCODE/geuvadis_drimseq_0_3_3_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_comparison_plots.Rout


##############################
### Positive controls for FDR = 0.1
##############################

### Prepare lists of PCR validated sQTLs and GWAS SNPs from GLiMMPS paper
### Run the "prepare_validated_sqtls.R" script




### Tables with adj p-values for validated sQTLs; Plots of this tables; Plots of expression of validated sQTLs

# GLIMMPS: PCR validated sQTLs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_fdr010' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_fdr010' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.1" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr.Rout



# GLIMMPS: GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_fdr010' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_fdr010' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.1" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas.Rout



# GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_fdr010' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_fdr010' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.1" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps.Rout



# GLIMMPS: sQTLs detected by GLIMMPS (~140)

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_glimmps.txt' plot_proportions=FALSE plot_tables=FALSE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_fdr010' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_fdr010' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.1" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps.Rout



# GEUVADIS: trQTLs detected within GEUVADIS project for the EUR population

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/geuvadis/geuvadis_valid_geuvadis.txt' plot_proportions=FALSE plot_tables=FALSE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_fdr010' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_fdr010' sqtlseeker_results='sqtlseeker_2_1_analysis' FDR=0.1" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis.Rout






### Plots the structure of the validated genes

# GLIMMPS: PCR validated sQTLs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_fdr010'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_summary.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_pcr.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_pcr.Rout

# GLIMMPS: GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt'  positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_fdr010'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_summary.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_gwas.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_gwas.Rout


# GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_fdr010'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_summary.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_gwas_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_gwas_glimmps.Rout




###################################################
## Comparison for FDR = 0.05 and sqtlseeker_2_1_analysis_drimseq_counts
###################################################

###### Comparisons for run with permutations: p-values from all the genes v2 - permutated p-values are saved under results$pvalue_perm



### Colors 
R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts'" $RCODE/colors.R $ROUT/colors.Rout


### Venn diagrams and upset plots with iCOBRA; histograms of p-values and nr. features 

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts' sqtlseeker_results='sqtlseeker_2_1_analysis_drimseq_counts' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_comparison_permutations.R $ROUT/geuvadis_drimseq_0_3_3_comparison_permutations.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_comparison_permutations.Rout


### For detected sQTLs: Distance to the closest exon, % within exons, mean gene expression, nr transcripts

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts' FDR=0.05 workers=10" $RCODE/geuvadis_drimseq_0_3_3_biol_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_biol_comparison_plots.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_biol_comparison_plots.Rout


### Plots of the overlap versus number of top ranked genes
### CAT (concordance at top) plots

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts' Overlaps_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotOverlaps.R' CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_comparison_plots.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_comparison_plots.Rout



##############################
### Positive controls for FDR = 0.05 and sqtlseeker_2_1_analysis_drimseq_counts
##############################

### Prepare lists of PCR validated sQTLs and GWAS SNPs from GLiMMPS paper
### Run the "prepare_validated_sqtls.R" script


### Tables with adj p-values for validated sQTLs; Plots of this tables; Plots of expression of validated sQTLs

# GLIMMPS: PCR validated sQTLs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_drimseq_counts' sqtlseeker_results='sqtlseeker_2_1_analysis_drimseq_counts' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr.Rout



# GLIMMPS: GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_drimseq_counts' sqtlseeker_results='sqtlseeker_2_1_analysis_drimseq_counts' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas.Rout



# GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt' plot_proportions=TRUE plot_tables=TRUE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_drimseq_counts' sqtlseeker_results='sqtlseeker_2_1_analysis_drimseq_counts' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps.Rout



# GLIMMPS: sQTLs detected by GLIMMPS (~140)

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_glimmps.txt' plot_proportions=FALSE plot_tables=FALSE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_drimseq_counts' sqtlseeker_results='sqtlseeker_2_1_analysis_drimseq_counts' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps.Rout



# GEUVADIS: trQTLs detected within GEUVADIS project for the EUR population

R32loc CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/geuvadis/geuvadis_valid_geuvadis.txt' plot_proportions=FALSE plot_tables=FALSE method_out='drimseq_0_3_3_analysis_permutations_all_genes' comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts' positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes_drimseq_counts' sqtlseeker_results='sqtlseeker_2_1_analysis_drimseq_counts' FDR=0.05" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis.Rout

tail $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis.Rout






##############################################################################

## Analysis on random BLOCK/SNP-gene pairs

##############################################################################

###################################################
## Run DRIMSeq analysis on random BLOCK/SNP-gene pairs 
###################################################


for n in {22..1}
do 

echo "${n}"

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 population='CEU' chr=${n}" $RCODE/geuvadis_drimseq_0_3_3_random_blocks.R $ROUT/geuvadis_drimseq_0_3_3_random_blocks_CEU_chr${n}.Rout

done

###################################################
## Run sQTLSeekeR analysis on random BLOCK/SNP-gene pairs
###################################################


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=5" $RCODE/geuvadis_sqtlseeker_2_1_random_blocks.R $ROUT/geuvadis_sqtlseeker_2_1_random_blocks.Rout



















