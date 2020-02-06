#===================================================================
#   File: Index_cSNP.sh
#   Directory code: /mnt/research/NMDL/2019_WB_MFM/cSNP  
#   Date: May 20, 2019
#   Description: Create index of cSNP for every 205e3 SNP inorder to run 
#           the Filtered_rst.R script in parallel.
#   Run: Rscript Index_cSNP.sh
#-------------------------------------------------------------------
#   Input files in directory:
#       /mnt/research/NMDL/2019_WB_MFM/cSNP/Function_Filter_cSNP
#
#   Input files:
#      VCF_MFM_WB.Rdata
#
#   Output file to directory:
#       /mnt/research/NMDL/2019_WB_MFM/cSNP
#
#   Output file:
#       Index_cSNP.txt
#===================================================================

#' Load required libraries
library(vcfR)

#' Clear Environment
rm(list=ls())

#' ### Load Data

#' VCF
vcf.Dir <- "/mnt/research/NMDL/2019_WB_MFM/cSNP/Function_Filter_cSNP"
load(paste(vcf.Dir, "VCF_MFM_WB.Rdata", sep="/"))

#' Number of cSNP
snp <- vcf@fix[,"ID"]
num <- length(snp)

#' Create a sequence to index file
x <- seq(from=1, to=num, by=205e3)
y <- c(x[-1] - 1, num)

#' Index
Ind <- data.frame(x,y)

#' Save index to text file
write.table(Ind, file=paste(getwd(), "Index_cSNP.txt", sep="/"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

