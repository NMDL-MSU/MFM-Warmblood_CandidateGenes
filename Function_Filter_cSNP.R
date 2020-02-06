
#' ### Description  
#' This scripts outlines the VCF file and the functions used to annotate and filter SNP. 
#'  
#' ***  
#' **Code:**  
#' Parent Directory:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP    
#'  
#' File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;Function_Filter_cSNP.R  
#'  
#' **Input files:**  
#' Directory/File  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/cSNP.vcf  
#'  
#' **Output files:**  
#'  
#' Directory:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/Function_Filter_cSNP    
#'  
#' Files:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;VCF_MFM_WB.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Functions_Annotate_SNP.Rdata  
#'  
#' Render R Script  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;Function_Filter_cSNP.qsub  
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

#' VCF cSNP File
vcf.Dir <- "/mnt/research/NMDL/2019_WB_MFM/cSNP"
vcf <- read.vcfR(paste(vcf.Dir, "cSNP.vcf", sep="/"), verbose = FALSE)

#' Add ID column to VCF data by concating the chromosome and position of each cSNP
vcf <- addID(vcf, sep = "_")
head(vcf)

#' Animal Ids
anim <- unlist(lapply(strsplit(colnames(vcf@gt)[-1], "_"), function(x) x[1]))
mfm <- c("9056", "9066", "12711", "12130", "12429", "9054", "12124", "12065")
anim <- data.frame(Anim=anim, MFM=ifelse(anim %in% mfm, "Afected", "Control"))
table(anim$MFM)


#' ### Table of VCF parameters  

#' #### INFO  
# Parameter
idx <- unlist(lapply(queryMETA(vcf, element="INFO"), function(x) x[1]))
param <- unlist(lapply(strsplit(idx, "="), function(x) x[3]))

# Description
idx <- unlist(lapply(queryMETA(vcf, element="INFO"), function(x) x[4]))
Desc <- unlist(lapply(strsplit(idx, "="), function(x) x[2]))

# Table
kable(data.frame(Parameter=param, Description=Desc), caption="VCF INFO")


#' #### FORMAT  
# Parameter  
idx <- unlist(lapply(queryMETA(vcf, element="FORMAT"), function(x) x[1]))
param <- unlist(lapply(strsplit(idx, "="), function(x) x[3]))

# Description
idx <- unlist(lapply(queryMETA(vcf, element="FORMAT"), function(x) x[4]))
Desc <- unlist(lapply(strsplit(idx, "="), function(x) x[2]))

# Table
kable(data.frame(Parameter=param, Description=Desc), caption="VCF INFO")


#' ### Filtering Function

#' Function to merge annotation (utility function called by filterSNP)
    annotate <- function(i = i, Info=Info, Ann=Ann, tmp=tmp, Des=Des, idx=idx, param=param){
        T <- strsplit(tmp[[i]], ";")[[1]]
        id <- unlist(lapply(strsplit(T, "="), function(x) x[1]))
        info <- do.call(cbind, lapply(strsplit(T, "="), function(x) x[2]))
        colnames(info) <- id

        if("INDEL" %in% id){
            info[,"INDEL"] <- "INDEL"
        }

        K <- data.frame(matrix(NA, nrow=1, ncol=17, dimnames=list(i, param[1:17])))
        for (j in param[1:17]){
            if (j %in% id){
                K[,j] <- info[,j]
            }
        }


        if("ANN" %in% id){
            ann <- strsplit(info[,"ANN"], "[|]")[[1]]
            ann <- matrix(ann, nrow=length(ann)/length(idx), ncol=length(idx), byrow=TRUE, 
                dimnames=list(NULL,idx))
            ann <- data.frame(SNP=rep(i, nrow(ann)), ann)
            ann$Allele <- gsub(",", "", ann$Allele)
        } else {
            ann <- matrix(NA, nrow=1, ncol=length(idx), dimnames=list(NULL, idx))
        }

        if("LOF" %in% id){
            lof <- gsub("[)]", "", info[,"LOF"])
            lof <- data.frame(LOF_Numb_transcripts_in_gene=strsplit(lof, "[|]")[[1]][3], 
                LOF_Percent_transcripts_affected=strsplit(lof, "[|]")[[1]][4])
        } else {
            lof <- data.frame(LOF_Numb_transcripts_in_gene=NA, LOF_Percent_transcripts_affected=NA)
        }

        if("NMD" %in% id){
            nmd <- gsub("[)]", "", info[,"NMD"])
            nmd <- data.frame(NMD_Numb_transcripts_in_gene=strsplit(nmd, "[|]")[[1]][3], 
                NMD_Percent_transcripts_affected=strsplit(nmd, "[|]")[[1]][4])
        } else {
            nmd <- data.frame(NMD_Numb_transcripts_in_gene=NA, NMD_Percent_transcripts_affected=NA)
        }

        Info[[i]] <- data.frame(data.frame(t(Des[i,])), K, lof, nmd)
        Ann[[i]] <- ann

        rst <- list(Info=Info[[i]], Ann=Ann[[i]])
        return(rst)
    }


#' The following function filters SNPs based on:

#' > Low quality: PASS filter (automatic taken from vcf)  
#' > `InDel`: removes them from downstream analysis (default=FALSE)  
#' > `MQSB`: Mann-Whitney U test of Mapping Quality vs Strand Bias p-value cutoff (default=0.05)  
#' > `Miss` : Genotype calls, filters SNP with less than `Miss` missing genotypes (default=2)  
#' > `TD` : Total allelic density: filters  SNP with less than `TD` reads (default=10)  

#' This function can be run in parallel using `cl` to specify the number of cores (default=20)

#' The required input is the vcf object of class vcfR (`vcf`)and a vector of SNP to filter 
#' (`SNP`, optional default is NULL). When specifying SNP make sure they match the ID in the vcf.
#' Default for SNP is to take the IDs directly from the vcf object (Note: computationally intensive) 

#'  The results are formated as list with four matrices and two lists:

#' > `Info`: SNP informtion (see INFO and FORMAT tables above for description)  
#' > `Ann` : A list containing SNP annotation (see INFO and FORMAT tables above for description)  
#' > `Geno`: Genotypes  
#' > `Ref`: Reference allele density  
#' > `Alt`: Alternative allele density  
#' > `FilteredInfo` : A list contating a summary table with the following information:  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`FiltSum` : Table with the number of SNP filtered out per category  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`FiltInfo` : Filteres SNP information  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`FiltAnn` : Filteres SNP annotation  


#' Filter variants Function
filterSNP <- function(vcf, SNP=NULL, InDel=FALSE, MQSB=0.05, Miss=2, TD=10, cl=20){
# Parameter
    idx <- unlist(lapply(queryMETA(vcf, element="INFO"), function(x) x[1]))
    param <- unlist(lapply(strsplit(idx, "="), function(x) x[3]))

# Description
    idx <- unlist(lapply(queryMETA(vcf, element="INFO"), function(x) x[4]))
    desc <- unlist(lapply(strsplit(idx, "="), function(x) x[2]))
    names(desc) <- param

# SNP and animal IDs
    if(is.null(SNP)){
        SNP <- vcf@fix[,"ID"]
    }
    ids <- vcf@fix[,"ID"]
    anim <- unlist(lapply(strsplit(colnames(vcf@gt)[-1], "_"), function(x) x[1]))

# SNP Information (chr, pos, etc)
    Des <-  vcf@fix[,c(1,2,4,5,6,7)]
    rownames(Des) <-  ids

# SNP Format
    tmp <- vcf@fix[,8]
    names(tmp) <- ids

# Annotation elements
    idx <- strsplit(desc[["ANN"]], "[|]")[[1]][-16]
    idx <- gsub("Functional annotations: '", "", idx)
    idx <- gsub(" ", "", idx)

# Merge SNP information and annotation
    Info <- list()
    Ann <- list()
    registerDoParallel(cl)
    rst <- foreach(i=SNP) %dopar% 
        annotate(i = i, Info=Info, Ann=Ann, tmp=tmp, Des=Des, idx=idx, param=param) 

# SNP Annotation
    INFO <- do.call(rbind, lapply(rst, function(x) x[["Info"]]))
    ANN <-  lapply(rst, function(x) x[["Ann"]])
    names(ANN) <- SNP

# Filter out low quality
    lowq <- SNP[!INFO[,"FILTER"] == "PASS"]

# Filter out INDELs
    if (InDel == TRUE){
        indel <- SNP[grep("INDEL", INFO[,"INDEL"])] 
    } else {
        indel <- NULL
    }

# Filter out SNP with significant mapping quality vs strand bias
    N <- as.numeric(INFO[,"MQSB"])
    names(N) <- rownames(INFO)
    N[is.na(N)] <- 0.04
    mqsb <- names(N)[N < MQSB]

# Filter based on genotype calls: less than animals with missing genotype
    gt <- vcf@gt
    rownames(gt) <- names(tmp)
    colnames(gt) <- c("FORMAT", anim)
    gt <- gt[rownames(INFO),]
    geno <- apply(gt, 2, function(x) unlist(lapply(strsplit(x, ":"), function(y) y[1])))
    geno <- geno[,-1]
    calls <- rownames(geno)[apply(geno, 1, function(x) length(grep("[.]", x))) > Miss]

# Filter based on base density: less than 10
    AD <- apply(gt, 2, function(x) unlist(lapply(strsplit(x, ":"), function(y) y[7])))
    AD <- AD[,-1]

# Alleleic Depth Reference
    ADR <- do.call(rbind, lapply(1:nrow(AD), 
        function(x) as.numeric(unlist(lapply(strsplit(AD[x,], ","), function(y) y[1])))))
    rownames(ADR) <- rownames(AD)
    colnames(ADR) <- colnames(AD)

# Allelic Depth Alternative
    ADA <- do.call(rbind, lapply(1:nrow(AD), 
        function(x) as.numeric(unlist(lapply(strsplit(AD[x,], ","), function(y) y[2])))))
    rownames(ADA) <- rownames(AD)
    colnames(ADA) <- colnames(AD)

# Total Allelic Density
    TAD <- ADR + ADA
    tad <- rownames(TAD)[apply(TAD, 1, function(x) sum(x < TD)) > Miss]

# Filtered SNP
    rmSNP <- unique(c(lowq, indel, mqsb, calls, tad))
    SNP <- SNP[!SNP %in% rmSNP]

# Table of filtered SNP
    Filt <- data.frame(LowQuality=length(lowq), InDel=length(indel), MQSB=length(mqsb), 
        GenoCall=length(calls), TAD=length(tad), TotalFiltered=length(rmSNP), 
        TotalRetained=length(SNP))

# Filtered SNP information
    FiltInfo <- INFO[rmSNP,]
    FiltAnn <- ANN[rmSNP]
    FilteredInfo <- list(FiltSum=Filt, FiltInfo=FiltInfo, FiltAnn=FiltAnn)

# Merge list of genotypes and allelic density
    rst <- list(Info=INFO[SNP,], Ann=ANN[SNP], Geno=geno[SNP,], Ref=ADR[SNP,], 
        Alt=ADA[SNP,], FilteredInfo=FilteredInfo)

    return(rst)
}

#' Save VCF object to file
save(vcf, file=paste(getwd(), "VCF_MFM_WB.Rdata", sep="/"))

#' Save functions
save(vcf, annotate, filterSNP, file=paste(getwd(), "Functions_Annotate_SNP.Rdata", sep="/"))

#' Session Information
sessionInfo()

#' ### Run R Script
#+ eval = FALSE
# htmlRunR
# Function_Filter_cSNP.R nodes=1,cpus-per-task=1,time=01:00:00,mem=10G +Function to filter and annotate SNP

