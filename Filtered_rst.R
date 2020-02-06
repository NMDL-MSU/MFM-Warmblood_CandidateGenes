#' ### Description  
#' Given the called cSNP from the MFM Warmblood RNA-seq data, filter out low  
#' quality cSNP and retain cSNP information, cSNP annotation and genotypes per horse
#' for index called SNP. 
#'  
#' ***  
#' **Code:**  
#' Parent Directory:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP    
#'  
#' File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;Filtered_rst.R  
#'  
#' **Input files:**  
#' Directory/File  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/Function_Filter_cSNP/VCF_MFM_WB.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/Function_Filter_cSNP/Functions_Annotate_SNP.Rdata  
#'  
#' **Output files:**  
#'  
#' Directory:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/Filtered_rst    
#'  
#' Files:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Filtered_rst.Rdata  
#'  
#' Render R Script  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;Filtered_rst.qsub  
#'  
#' ***  

#' ### R Environment

#' Load required libraries
library(vcfR)
library(foreach)
library(doParallel)
library(knitr)

#' Clear Environment
rm(list=ls())


#' ### Load Data

#' VCF
vcf.Dir <- "/mnt/research/NMDL/2019_WB_MFM/cSNP/Function_Filter_cSNP"
load(paste(vcf.Dir, "VCF_MFM_WB.Rdata", sep="/"))

#' Functions
load(paste(vcf.Dir, "Functions_Annotate_SNP.Rdata", sep="/"))

#' Parameters to filter variants:
ind <- index
cSNP <- vcf@fix[,"ID"]
InDel <- TRUE
MQSB <- 0.05 
Miss <- 2
TD <- 10
cl <- 20

#' Run filter and annotate function
system.time({
    rst <- filterSNP(vcf, SNP=cSNP[ind], InDel=TRUE)
})

#' Save results to file
save(rst, file=paste(getwd(), "Filtered_rst.Rdata", sep="/"))

#' Session Information
sessionInfo()

#' ### Run R Script
#+ eval = FALSE
# htmlRunR
# Filtered_rst.R nodes=1,cpus-per-task=21,time=04:00:00,mem=300G +WB MFM index cSNP Filter

