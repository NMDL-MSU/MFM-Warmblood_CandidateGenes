---
title: Function to filter and annotate SNP
author: Deborah Velez-Irizarry
date: Fri May 17 14:29:56 EDT 2019
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description  
This scripts outlines the VCF file and the functions used to annotate and filter SNP. 
 
***  
**Code:**  
Parent Directory:  
 
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP    
 
File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;Function_Filter_cSNP.R  
 
**Input files:**  
Directory/File  
 
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/cSNP.vcf  
 
**Output files:**  
 
Directory:  
 
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/Function_Filter_cSNP    
 
Files:  
 
>&nbsp;&nbsp;&nbsp;&nbsp;VCF_MFM_WB.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;Functions_Annotate_SNP.Rdata  
 
Render R Script  
 
> &nbsp;&nbsp;&nbsp;&nbsp;Function_Filter_cSNP.qsub  
 
***  
### R Environment
Load required libraries


```r
library(vcfR)
```

```
## 
##    *****       ***   vcfR   ***       *****
##    This is vcfR 1.8.0 
##      browseVignettes('vcfR') # Documentation
##      citation('vcfR') # Citation
##    *****       *****      *****       *****
```

```r
library(foreach)
library(doParallel)
```

```
## Loading required package: iterators
```

```
## Loading required package: parallel
```

```r
library(knitr)
```

Clear Environment


```r
rm(list=ls())
```

### Load Data
VCF cSNP File


```r
vcf.Dir <- "/mnt/research/NMDL/2019_WB_MFM/cSNP"
vcf <- read.vcfR(paste(vcf.Dir, "cSNP.vcf", sep="/"), verbose = FALSE)
```

Add ID column to VCF data by concating the chromosome and position of each cSNP


```r
vcf <- addID(vcf, sep = "_")
head(vcf)
```

```
## [1] "***** Object of class 'vcfR' *****"
## [1] "***** Meta section *****"
## [1] "##fileformat=VCFv4.2"
## [1] "##FILTER=<ID=PASS,Description=\"All filters passed\">"
## [1] "##bcftoolsVersion=1.9-64-g28bcc56+htslib-1.9-52-g6e86e38"
## [1] "##bcftoolsCommand=mpileup -Ou -C50 -E -Q25 -a DV,AD,ADF,ADR,SP -f /m [Truncated]"
## [1] "##reference=file:///mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughb [Truncated]"
## [1] "##contig=<ID=chr1,length=188260577>"
## [1] "First 6 rows."
## [1] 
## [1] "***** Fixed section *****"
##      CHROM  POS   ID           REF ALT QUAL      FILTER   
## [1,] "chr1" "229" "chr1_229_1" "C" "G" "513.0"   "PASS"   
## [2,] "chr1" "233" "chr1_233_2" "G" "T" "162.0"   "PASS"   
## [3,] "chr1" "265" "chr1_265_3" "T" "G" "763.0"   "PASS"   
## [4,] "chr1" "311" "chr1_311_4" "T" "C" "758.0"   "PASS"   
## [5,] "chr1" "363" "chr1_363_5" "G" "A" "999.0"   "PASS"   
## [6,] "chr1" "377" "chr1_377_6" "T" "C" "14.9975" "LowQual"
## [1] 
## [1] "***** Genotype section *****"
##      FORMAT                   12065_uniq_sorted.bam         
## [1,] "GT:PL:DV:SP:ADF:ADR:AD" "0/1:23,0,167:1:0:1,1:4,0:5,1"
## [2,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,18,208:0:0:2,0:4,0:6,0"
## [3,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,3,41:0:0:0,0:1,0:1,0"  
## [4,] "GT:PL:DV:SP:ADF:ADR:AD" "0/1:33,0,31:1:0:0,0:1,1:1,1" 
## [5,] "GT:PL:DV:SP:ADF:ADR:AD" "0/1:62,0,79:2:0:1,0:1,2:2,2" 
## [6,] "GT:PL:DV:SP:ADF:ADR:AD" "0/0:0,12,150:0:0:1,0:3,0:4,0"
##      12124_uniq_sorted.bam          12130_uniq_sorted.bam         
## [1,] "0/0:0,12,150:0:0:1,0:3,0:4,0" "1/1:165,15,0:5:0:0,3:0,2:0,5"
## [2,] "0/1:26,0,113:1:0:1,0:2,1:3,1" "0/0:0,15,170:0:0:3,0:2,0:5,0"
## [3,] "0/1:68,0,32:2:0:0,1:1,1:1,2"  "1/1:129,12,0:4:0:0,3:0,1:0,4"
## [4,] "0/1:24,0,110:1:0:4,1:0,0:4,1" "1/1:39,3,0:1:0:0,0:0,1:0,1"  
## [5,] "0/0:0,27,251:0:0:7,0:2,0:9,0" "1/1:230,24,0:8:0:0,4:0,4:0,8"
## [6,] "0/0:0,21,222:0:0:5,0:2,0:7,0" "0/0:0,27,255:0:0:5,0:4,0:9,0"
##      12212_uniq_sorted.bam           12216_uniq_sorted.bam           
## [1,] "0/1:145,0,192:5:0:3,2:3,3:6,5" "0/1:98,0,133:3:0:1,1:3,2:4,3"  
## [2,] "0/0:0,33,255:0:0:6,0:5,0:11,0" "0/0:0,21,223:0:0:2,0:5,0:7,0"  
## [3,] "0/1:92,0,99:3:0:1,2:2,1:3,3"   "0/1:54,0,110:2:4:2,0:1,2:3,2"  
## [4,] "0/1:81,0,183:3:0:2,1:4,2:6,3"  "0/1:64,0,70:2:0:1,1:1,1:2,2"   
## [5,] "0/1:156,0,210:6:6:3,5:4,1:7,6" "0/1:212,0,142:8:0:4,5:1,3:5,8" 
## [6,] "0/0:0,39,255:0:0:8,0:5,0:13,0" "0/0:0,42,255:0:0:10,0:4,0:14,0"
## [1] "First 6 columns only."
## [1] 
## [1] "Unique GT formats:"
## [1] "GT:PL:DV:SP:ADF:ADR:AD"
## [1]
```

Animal Ids


```r
anim <- unlist(lapply(strsplit(colnames(vcf@gt)[-1], "_"), function(x) x[1]))
mfm <- c("9056", "9066", "12711", "12130", "12429", "9054", "12124", "12065")
anim <- data.frame(Anim=anim, MFM=ifelse(anim %in% mfm, "Afected", "Control"))
table(anim$MFM)
```

```
## 
## Afected Control 
##       8       8
```

### Table of VCF parameters  
#### INFO  


```r
# Parameter
idx <- unlist(lapply(queryMETA(vcf, element="INFO"), function(x) x[1]))
param <- unlist(lapply(strsplit(idx, "="), function(x) x[3]))

# Description
idx <- unlist(lapply(queryMETA(vcf, element="INFO"), function(x) x[4]))
Desc <- unlist(lapply(strsplit(idx, "="), function(x) x[2]))

# Table
kable(data.frame(Parameter=param, Description=Desc), caption="VCF INFO")
```



|Parameter |Description                                                                                                                                                                                                                                                                                                                                             |
|:---------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|INDEL     |Indicates that the variant is an INDEL.                                                                                                                                                                                                                                                                                                                 |
|IDV       |Maximum number of reads supporting an indel                                                                                                                                                                                                                                                                                                             |
|IMF       |Maximum fraction of reads supporting an indel                                                                                                                                                                                                                                                                                                           |
|DP        |Raw read depth                                                                                                                                                                                                                                                                                                                                          |
|VDB       |Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)"                                                                                                                                                                                                                                                           |
|RPB       |Mann-Whitney U test of Read Position Bias (bigger is better)                                                                                                                                                                                                                                                                                            |
|MQB       |Mann-Whitney U test of Mapping Quality Bias (bigger is better)                                                                                                                                                                                                                                                                                          |
|BQB       |Mann-Whitney U test of Base Quality Bias (bigger is better)                                                                                                                                                                                                                                                                                             |
|MQSB      |Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)                                                                                                                                                                                                                                                                                |
|SGB       |Segregation based metric.                                                                                                                                                                                                                                                                                                                               |
|MQ0F      |Fraction of MQ0 reads (smaller is better)                                                                                                                                                                                                                                                                                                               |
|ICB       |Inbreeding Coefficient Binomial test (bigger is better)                                                                                                                                                                                                                                                                                                 |
|HOB       |Bias in the number of HOMs number (smaller is better)                                                                                                                                                                                                                                                                                                   |
|AC        |Allele count in genotypes for each ALT allele                                                                                                                                                                                                                                                                                                           |
|AN        |Total number of alleles in called genotypes                                                                                                                                                                                                                                                                                                             |
|DP4       |Number of high-quality ref-forward                                                                                                                                                                                                                                                                                                                      |
|MQ        |Average mapping quality                                                                                                                                                                                                                                                                                                                                 |
|ANN       |Functional annotations: 'Allele &#124; Annotation &#124; Annotation_Impact &#124; Gene_Name &#124; Gene_ID &#124; Feature_Type &#124; Feature_ID &#124; Transcript_BioType &#124; Rank &#124; HGVS.c &#124; HGVS.p &#124; cDNA.pos / cDNA.length &#124; CDS.pos / CDS.length &#124; AA.pos / AA.length &#124; Distance &#124; ERRORS / WARNINGS / INFO' |
|LOF       |Predicted loss of function effects for this variant. Format: 'Gene_Name &#124; Gene_ID &#124; Number_of_transcripts_in_gene &#124; Percent_of_transcripts_affected'                                                                                                                                                                                     |
|NMD       |Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name &#124; Gene_ID &#124; Number_of_transcripts_in_gene &#124; Percent_of_transcripts_affected'                                                                                                                                                                              |

#### FORMAT  


```r
# Parameter  
idx <- unlist(lapply(queryMETA(vcf, element="FORMAT"), function(x) x[1]))
param <- unlist(lapply(strsplit(idx, "="), function(x) x[3]))

# Description
idx <- unlist(lapply(queryMETA(vcf, element="FORMAT"), function(x) x[4]))
Desc <- unlist(lapply(strsplit(idx, "="), function(x) x[2]))

# Table
kable(data.frame(Parameter=param, Description=Desc), caption="VCF INFO")
```



|Parameter |Description                                |
|:---------|:------------------------------------------|
|PL        |List of Phred-scaled genotype likelihoods  |
|DV        |Number of high-quality non-reference bases |
|SP        |Phred-scaled strand bias P-value           |
|AD        |Allelic depths                             |
|ADF       |Allelic depths on the forward strand       |
|ADR       |Allelic depths on the reverse strand       |
|GT        |Genotype                                   |

### Filtering Function
Function to merge annotation (utility function called by filterSNP)


```r
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
```

The following function filters SNPs based on:
> Low quality: PASS filter (automatic taken from vcf)  
> `InDel`: removes them from downstream analysis (default=FALSE)  
> `MQSB`: Mann-Whitney U test of Mapping Quality vs Strand Bias p-value cutoff (default=0.05)  
> `Miss` : Genotype calls, filters SNP with less than `Miss` missing genotypes (default=2)  
> `TD` : Total allelic density: filters  SNP with less than `TD` reads (default=10)  
This function can be run in parallel using `cl` to specify the number of cores (default=20)
The required input is the vcf object of class vcfR (`vcf`)and a vector of SNP to filter 
(`SNP`, optional default is NULL). When specifying SNP make sure they match the ID in the vcf.
Default for SNP is to take the IDs directly from the vcf object (Note: computationally intensive) 
 The results are formated as list with four matrices and two lists:
> `Info`: SNP informtion (see INFO and FORMAT tables above for description)  
> `Ann` : A list containing SNP annotation (see INFO and FORMAT tables above for description)  
> `Geno`: Genotypes  
> `Ref`: Reference allele density  
> `Alt`: Alternative allele density  
> `FilteredInfo` : A list contating a summary table with the following information:  
> &nbsp;&nbsp;&nbsp;&nbsp;`FiltSum` : Table with the number of SNP filtered out per category  
> &nbsp;&nbsp;&nbsp;&nbsp;`FiltInfo` : Filteres SNP information  
> &nbsp;&nbsp;&nbsp;&nbsp;`FiltAnn` : Filteres SNP annotation  
Filter variants Function


```r
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
```

Save VCF object to file


```r
save(vcf, file=paste(getwd(), "VCF_MFM_WB.Rdata", sep="/"))
```

Save functions


```r
save(vcf, annotate, filterSNP, file=paste(getwd(), "Functions_Annotate_SNP.Rdata", sep="/"))
```

Session Information


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /opt/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas_sandybridgep-r0.3.1.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] doParallel_1.0.11 iterators_1.0.10  foreach_1.4.4     vcfR_1.8.0       
## [5] knitr_1.20       
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.1        cluster_2.0.7-1   magrittr_1.5     
##  [4] MASS_7.3-51.1     viridisLite_0.3.0 ape_5.1          
##  [7] lattice_0.20-38   highr_0.7         stringr_1.4.0    
## [10] tools_3.5.1       grid_3.5.1        nlme_3.1-137     
## [13] mgcv_1.8-24       vegan_2.5-4       permute_0.9-5    
## [16] Matrix_1.2-14     codetools_0.2-15  evaluate_0.12    
## [19] pinfsc50_1.1.0    stringi_1.4.3     compiler_3.5.1   
## [22] memuse_4.0-0
```

### Run R Script


```r
# htmlRunR
# Function_Filter_cSNP.R nodes=1,cpus-per-task=1,time=01:00:00,mem=10G +Function to filter and annotate SNP
```

