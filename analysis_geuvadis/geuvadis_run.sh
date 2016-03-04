#!/bin/bash
## Define paths to software and reference files

RCODE=/home/gosia/R/drimseq_paper/analysis_geuvadis
RWD=/home/Shared/data/seq/geuvadis
ROUT=$RWD/Rout

mkdir $ROUT

###################################################
## Run "prepare_data.R"
###################################################


##############################
### Colors
##############################

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD/drimseq_0_3_3_comparison'" $RCODE/colors.R $ROUT/colors.Rout


###################################################
## Run DRIMSeq analysis
###################################################


for n in {22..1}
do 

echo "${n}"

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4 population='CEU' chr=${n}" $RCODE/geuvadis_drimseq_0_3_3_run.R $ROUT/geuvadis_drimseq_0_3_3_run_CEU_chr${n}.Rout

done


###################################################
## Run sQTLSeekeR analysis
###################################################


R31 CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=4" $RCODE/geuvadis_sqtlseeker_2_1_run.R $ROUT/geuvadis_sqtlseeker_2_1_run.Rout


###################################################
## Comparison 
###################################################

### Venn diagrams 

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU'" $RCODE/geuvadis_drimseq_0_3_3_comparison.R $ROUT/geuvadis_drimseq_0_3_3_comparison_run.Rout


### Plots of the overlap versus number of top ranked genes

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' Overlaps_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotOverlaps.R'" $RCODE/geuvadis_drimseq_0_3_3_comparison_plots.R $ROUT/geuvadis_drimseq_0_3_3_comparison_plots.Rout


##############################
### Positive controls
##############################

### Prepare lists of PCR validated sQTLs and GWAS SNPs from GLiMMPS paper
### Run the "prepare_validated_sqtls.R" script


### Tables with adj p-values for validated sQTLs; Plots of expression of validated sQTLs

# GLIMMPS: PCR validated sQTLs

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt' plot_proportions=TRUE" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_pcr.Rout

# GLIMMPS: GWAS SNPs

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt' plot_proportions=TRUE" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas.Rout


# GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt' plot_proportions=TRUE" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_gwas_glimmps.Rout


# GLIMMPS: sQTLs detected by GLIMMPS (~140)

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/glimmps/glimmps_valid_glimmps.txt' plot_proportions=FALSE" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_glimmps.Rout


# GEUVADIS: trQTLs detected within GEUVADIS project for the EUR population

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' valid_path='data/validation/geuvadis/geuvadis_valid_geuvadis.txt' plot_proportions=FALSE" $RCODE/geuvadis_drimseq_0_3_3_positive_controls.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_geuvadis.Rout



### Plots of tables with adj p-values for validated sQTLs; Plots of gene structure

# GLIMMPS: PCR validated sQTLs

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_pcr.txt'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_summary.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_pcr.Rout

# GLIMMPS: GWAS SNPs

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_gwas.txt'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_summary.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_gwas.Rout


# GLIMMPS: GLIMMPS SNPs that are linked to GWAS SNPs

R32 CMD BATCH --no-save --no-restore "--args rwd='$RWD' population='CEU' path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf' valid_path='data/validation/glimmps/glimmps_valid_gwas_glimmps.txt'" $RCODE/geuvadis_drimseq_0_3_3_positive_controls_summary.R $ROUT/geuvadis_drimseq_0_3_3_positive_controls_summary_gwas_glimmps.Rout


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




##############################################################################

## Run DRIMSeq analysis using F test - R32dev

##############################################################################


for n in 22
do 

echo "${n}"

R32dev CMD BATCH --no-save --no-restore "--args rwd='$RWD' workers=6 population='CEU' chr=${n}" $RCODE/geuvadis_drimseq_0_3_3_run_f.R $ROUT/geuvadis_drimseq_0_3_3_run_f_CEU_chr${n}.Rout

done


















