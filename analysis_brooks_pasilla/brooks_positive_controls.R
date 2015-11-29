######################################################
## ----- brooks_positive_controls
## <<brooks_positive_controls.R>>

# BioC 3.1
# Created 29 Nov 2015 

##############################################################################

library(plyr)
library(ggplot2)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/brooks_pasilla'


##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)



##############################################################################

setwd(rwd)

comparison_out <- "drimseq_0_3_1_positive_controls/"
dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)

out_dir <- paste0(comparison_out)

##############################################################################

