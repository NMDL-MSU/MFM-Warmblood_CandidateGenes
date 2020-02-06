---
title: Genotypes Public Data
author: Deborah Velez-Irizarry
date: Mon Apr 1 16:21:02 EDT 2019
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---


```r
### Description:  
```

Extract genotypes of interest from public data sources.  
[Warmblood PRJEB23301](https://www.ebi.ac.uk/ena/data/view/PRJEB23301)
[Mix Breeds PRJEB28306](https://www.ebi.ac.uk/ena/data/view/PRJEB28306)
[Quarter Horses PRJEB30116](https://www.ebi.ac.uk/ena/data/view/PRJEB30116)
 
***  
**Code:**  
Parent Directory:  
 
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/Public_Data  
 
File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;GenotypesPublicData.R  
 
**Input files:**  
Directories  
 
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/Public_Data  
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/DNA/VariantCalling  
  
Files:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;379_RAO_HD_Warmbloods.ne1000.vcf  
>&nbsp;&nbsp;&nbsp;&nbsp;horses.88.vars.flt.pass.ebi.vcf  
>&nbsp;&nbsp;&nbsp;&nbsp;cSNP.vcf  
 
**Output files:**  
 
Directory:  
 
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/Public_Data/GenotypesPublicData  
 
Files:  
 
>&nbsp;&nbsp;&nbsp;&nbsp;TB.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;genoTB.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;dataTB.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;snp_data_TB.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;WB379.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;genoWB379.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;snp_data_WB379.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Mix88.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;genoMix88.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;snp_data_Mix88.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;QH.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;genoQH.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;snp_data_QH.txt  
 
Render R Script  
 
> &nbsp;&nbsp;&nbsp;&nbsp;GenotypesPublicData.qsub  
 
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

Clear Environment


```r
rm(list=ls())
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] vcfR_1.8.0 knitr_1.21
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.1        lattice_0.20-38   ape_5.1          
##  [4] viridisLite_0.3.0 permute_0.9-4     MASS_7.3-51.1    
##  [7] grid_3.5.1        nlme_3.1-137      magrittr_1.5     
## [10] evaluate_0.12     highr_0.7         stringi_1.2.3    
## [13] vegan_2.5-2       Matrix_1.2-14     pinfsc50_1.1.0   
## [16] tools_3.5.1       stringr_1.3.1     xfun_0.4         
## [19] parallel_3.5.1    compiler_3.5.1    cluster_2.0.7-1  
## [22] mgcv_1.8-24
```

### Load Data
VCF cSNP File


```r
Dir1 <- "/mnt/research/NMDL/2019_WB_MFM/Public_Data"
Dir2 <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/DNA/VariantCalling"
vcf.files <- c(paste(Dir1, list.files(Dir1)[grep("vcf", list.files(Dir1))], sep="/"),
    paste(Dir2, list.files(Dir2)[grep("cSNP.vcf", list.files(Dir2))], sep="/"))
names(vcf.files) <- c("WB379", "QH", "Mix88", "TB")
vcf.files
```

```
##                                                                              WB379 
##      "/mnt/research/NMDL/2019_WB_MFM/Public_Data/379_RAO_HD_Warmbloods.ne1000.vcf" 
##                                                                                 QH 
##                     "/mnt/research/NMDL/2019_WB_MFM/Public_Data/ecab_variants.vcf" 
##                                                                              Mix88 
##       "/mnt/research/NMDL/2019_WB_MFM/Public_Data/horses.88.vars.flt.pass.ebi.vcf" 
##                                                                                 TB 
## "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/DNA/VariantCalling/cSNP.vcf"
```

SNP of interest


```r
posSNP <- read.table("/mnt/research/NMDL/2019_WB_MFM/Public_Data/snp_of_interest.txt", header=TRUE)
```

### Read VCF files into R: Thoroughbred Horses


```r
TB <- read.vcfR(vcf.files["TB"], verbose = FALSE)
```

Extract VCF INFO


```r
info.snp <- data.frame(getFIX(TB))
info.snp$POS <- as.numeric(as.character(info.snp$POS))
nrow(info.snp)
```

```
## [1] 2355791
```

Extract SNPs within genes of interest


```r
SNP <- lapply(1:nrow(posSNP), function(x) 
    data.frame(gene=posSNP$gene[x], info.snp[as.character(info.snp$CHROM) == as.character(posSNP$chr[x]) & 
        info.snp$POS == posSNP$EquCab3[x],]))
SNP <- do.call(rbind, SNP)
SNP
```

```
##            gene CHROM       POS   ID REF ALT    QUAL  FILTER
## 12526      BAG3  chr1  12966053 <NA>   C   T     210    PASS
## 81585      LDB3  chr1  84589937 <NA>   T   C     999    PASS
## 469438     FLNC  chr4  83840299 <NA>   G   A     999    PASS
## 496124   DNAJB6  chr4 107802363 <NA>   G   T     999    PASS
## 544211     LMNA  chr5  38652136 <NA>   A   C     999    PASS
## 608220      DES  chr6   8696183 <NA>   T   G     210    PASS
## 648785  PYROXD1  chr6  48903284 <NA>   A   G     999    PASS
## 648874  PYROXD1  chr6  48918053 <NA>   G   A     999    PASS
## 648875  PYROXD1  chr6  48918068 <NA>   G   A     999    PASS
## 648935  PYROXD1  chr6  48924749 <NA>   G   C     913    PASS
## 957534     PLEC  chr9  84701301 <NA>   C   T     778    PASS
## 957535     PLEC  chr9  84701459 <NA>   C   T     210    PASS
## 957537     PLEC  chr9  84701931 <NA>   C   T     683    PASS
## 957541     PLEC  chr9  84702701 <NA>   G   A     447    PASS
## 957542     PLEC  chr9  84702893 <NA>   A   G     210    PASS
## 957547     PLEC  chr9  84703445 <NA>   C   T     999    PASS
## 957548     PLEC  chr9  84703610 <NA>   G   A     210    PASS
## 957554     PLEC  chr9  84704450 <NA>   T   C     999    PASS
## 957557     PLEC  chr9  84705316 <NA>   C   T     687    PASS
## 957590     PLEC  chr9  84715932 <NA>   C   A     447    PASS
## 1307578  SQSTM1 chr14   1899276 <NA>   G   A     999    PASS
## 1307654  SQSTM1 chr14   1910307 <NA>   G   A     210    PASS
## 1341857    MYOT chr14  37814680 <NA>   G   C     210    PASS
## 1341877    MYOT chr14  37818807 <NA>   A   G     999    PASS
## 1341878    MYOT chr14  37818823 <NA>   A   G     999    PASS
## 1548240      KY chr16  71665337 <NA>   C   T 27.6939 LowQual
```

Extract Genotypes


```r
geno <- extract.gt(TB, return.alleles=FALSE)
colnames(geno) <- sapply(strsplit(colnames(geno), "_"), function(x) x[[1]][1])
dim(geno)
```

```
## [1] 2355791      79
```

```r
genoTB <- geno[as.numeric(rownames(SNP)),] 
dim(genoTB)
```

```
## [1] 26 79
```

Animal Information RER Thoroughbred


```r
RER <- read.table("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq/RER_Thoroughbred_Animal_Information.txt", 
    header=TRUE, sep="\t", row.names=1)
RER$ID <- as.character(RER$ID)
rownames(RER) <- RER$ID
```

Animal Information for Glycogen Thoroughbred


```r
Gly <- read.table("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/Glycogen_Project_Information.txt",
    header=TRUE, sep="\t")

# Retain information on sequenced animals 
Gly <- Gly[!is.na(Gly$MSMS_Plate),]
Gly$TimePoint <- factor(as.character(Gly$TimePoint), exclude="Rep48h", 
    levels=c("Pre", "Depl", "Rep24h", "Rep72h"))
Gly$Diet <- factor(as.character(Gly$Diet), exclude="DietGoldenMax")
rownames(Gly) <- paste("G", Gly$MSMS_ID, sep="")
head(Gly)
```

```
##    Animal Period      Diet DietStarch Horse TimePoint DateTrial MSMS_Plate
## G1  A7934      1       Fat        Low    Pi       Pre 3/27/2012          1
## G2  A7934      1       Fat        Low    Pi      Depl 3/30/2012          1
## G3  A7934      1       Fat        Low    Pi    Rep24h 3/31/2012          1
## G4  A7934      1       Fat        Low    Pi    Rep72h  4/2/2012          1
## G5  A7932      1 SweetFeed       High  King       Pre 3/27/2012          1
## G6  A7932      1 SweetFeed       High  King      Depl 3/30/2012          1
##    MSMS_ID GlycoMN GlycoKA
## G1       1 128.855     130
## G2       2  96.115     138
## G3       3 104.645     176
## G4       4  58.921     112
## G5       5  95.094     108
## G6       6  89.002     142
```

```r
# Animal Ids for Glycogen Thoroughbreds
idx <- lapply(unique(Gly$Animal), function(x) rownames(Gly[Gly$Animal == x,]))
names(idx) <- unique(Gly$Animal)
idx
```

```
## $A7934
## [1] "G1"  "G2"  "G3"  "G4"  "G33" "G34" "G35" "G36"
## 
## $A7932
## [1] "G5"  "G6"  "G7"  "G8"  "G21" "G22" "G23" "G24"
## 
## $A7935
## [1] "G9"  "G10" "G11" "G12" "G37" "G38" "G39" "G40"
## 
## $A7933
## [1] "G13" "G14" "G15" "G16" "G25" "G26" "G27" "G28"
## 
## $A7937
## [1] "G17" "G18" "G19" "G20" "G29" "G30" "G31" "G32"
```

Check that genotypes are consistent across replicates


```r
ck <- lapply(idx, function(x) unlist(apply(genoTB[,x], 1, function(y) sum(!y %in% y[1]))))
lapply(ck, function(x) x[x>0])
```

```
## $A7934
## chr4_107802363_496124  chr9_84704450_957554 
##                     1                     1 
## 
## $A7932
## named integer(0)
## 
## $A7935
## chr4_107802363_496124  chr6_48903284_648785  chr6_48918068_648875 
##                     1                     1                     1 
##  chr9_84704450_957554 
##                     1 
## 
## $A7933
## chr16_71665337_1548240 
##                      1 
## 
## $A7937
## named integer(0)
```

Remove repeated animals from genotype matrix


```r
un <- unlist(lapply(idx, function(x) x[1]))
unG <-  genoTB[,(ncol(geno)-40):ncol(geno)][,un]
colnames(unG) <- names(un)
```

RER genotypes


```r
rerG <- genoTB[,colnames(genoTB) %in% RER$ID]
```

TB Genotype Matrix


```r
genoTB <- cbind(unG, rerG)
dim(genoTB)
```

```
## [1] 26 28
```

SNP Info for TB


```r
SNPinfo.TB <- SNP
```

Data matrix for TB


```r
dataTB <- rbind(data.frame(ID=Gly[un, "Animal"], Dx=rep("Control", length(un)), Breed=rep("Thoroughbred", length(un))),
    data.frame(RER[colnames(rerG),c("ID", "Dx")], Breed=rep("Thoroughbred", ncol(rerG))))
rownames(dataTB) <- NULL
```

Save genotypes for Thoroughbreds


```r
SNP.TB <- SNP
save(TB, genoTB, SNP.TB, dataTB, file="TB.Rdata")
write.table(genoTB, file="genotypesTB.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(SNP.TB, file="snp_data_TB.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(dataTB, file="dataTB.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
```

### Read VCF files into R: Warmblood Horses 379


```r
WB379 <- read.vcfR(vcf.files["WB379"], verbose = FALSE)
```

Extract VCF INFO


```r
info.snp <- data.frame(getFIX(WB379))
info.snp$POS <- as.numeric(as.character(info.snp$POS))
nrow(info.snp)
```

```
## [1] 1926709
```

Extract SNPs within genes of interest


```r
SNP <- lapply(1:nrow(posSNP), function(x) 
    data.frame(info.snp[as.character(info.snp$CHROM) == as.character(posSNP$chr[x]) & 
        info.snp$POS == posSNP$EquCab2[x],]))
idx <- unlist(lapply(SNP, nrow)) > 0
SNP <- data.frame(gene=posSNP$gene[idx], do.call(rbind, SNP))
SNP
```

```
##         gene CHROM      POS                   ID REF ALT QUAL FILTER
## 1831004 PLEC  chr9 82500455 MNEc.2.9.82500455.PC   C   T <NA>   PASS
## 1831006 PLEC  chr9 82501697 MNEc.2.9.82501697.PC   G   A <NA>   PASS
```

Extract Genotypes


```r
geno <- extract.gt(WB379, return.alleles=FALSE)
colnames(geno) <- sapply(strsplit(colnames(geno), "_"), function(x) x[[1]][1])
dim(geno)
```

```
## [1] 1926709     379
```

```r
genoWB379 <- geno[as.numeric(rownames(SNP)),] 
dim(genoWB379)
```

```
## [1]   2 379
```

Save genotypes for Warmblood


```r
SNP.WB379 <- SNP
save(WB379, genoWB379, SNP.WB379, file="WB379.Rdata")
write.table(genoWB379, file="genotypesWB379.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(SNP.WB379, file="snp_data_WB379.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
```

### Read VCF files into R: Mixed Breed Horses 88


```r
system.time(
    Mix88 <- read.vcfR(vcf.files["Mix88"], verbose = FALSE)
)
```

```
##     user   system  elapsed 
## 3184.129  310.792 4029.759
```

Extract VCF INFO


```r
info.snp <- data.frame(getFIX(Mix88))
info.snp$POS <- as.numeric(as.character(info.snp$POS))
nrow(info.snp)
```

```
## [1] 26128714
```

Extract SNPs within genes of interest


```r
SNP <- lapply(1:nrow(posSNP), function(x) 
    data.frame(info.snp[as.character(info.snp$CHROM) == as.character(posSNP$chr[x]) & 
        info.snp$POS == posSNP$EquCab3[x],]))
idx <- unlist(lapply(SNP, nrow)) > 0
SNP <- data.frame(gene=posSNP$gene[idx], do.call(rbind, SNP))
SNP
```

```
##             gene CHROM       POS   ID  REF    ALT     QUAL FILTER
## 7512658     BAG3  chr1  12966053 <NA> CGGG TGGG,C  1095.47   PASS
## 8168533     LDB3  chr1  84589937 <NA>    T      C 10548.73   PASS
## 18221518    FLNC  chr4  83840299 <NA>    G      A  1468.39   PASS
## 18469675  DNAJB6  chr4 107802363 <NA>    G      T  5154.83   PASS
## 18872296    LMNA  chr5  38652136 <NA>    A      C 12271.06   PASS
## 19889265 PYROXD1  chr6  48903284 <NA>    A      G 50141.19   PASS
## 19889475 PYROXD1  chr6  48918053 <NA>    G      A 19511.74   PASS
## 19889476 PYROXD1  chr6  48918068 <NA>    G      A  9880.86   PASS
## 19889559 PYROXD1  chr6  48924749 <NA>    G      C  2989.80   PASS
## 22968334    PLEC  chr9  84701301 <NA>    C      T  2275.08   PASS
## 22968336    PLEC  chr9  84701459 <NA>    C      T   559.98   PASS
## 22968343    PLEC  chr9  84701931 <NA>    C      T  1321.69   PASS
## 22968364    PLEC  chr9  84702701 <NA>    G      A  2897.64   PASS
## 22968365    PLEC  chr9  84702893 <NA>    A      G  2471.36   PASS
## 22968369    PLEC  chr9  84703445 <NA>    C      T 12875.47   PASS
## 22968416    PLEC  chr9  84704450 <NA>    T      C  7449.75   PASS
## 22968429    PLEC  chr9  84705316 <NA>    C      T  3366.76   PASS
## 22968585    PLEC  chr9  84715932 <NA>    C      A   898.91   PASS
## 2526992   SQSTM1 chr14   1899276 <NA>    G      A   561.32   PASS
## 2527132   SQSTM1 chr14   1910307 <NA>    G      A   950.80   PASS
## 2856860     MYOT chr14  37814680 <NA>    G      C   249.21   PASS
## 2856890     MYOT chr14  37818807 <NA>    A      G  9652.82   PASS
## 2856891     MYOT chr14  37818823 <NA>    A      G  7417.87   PASS
```

Extract Genotypes


```r
geno <- extract.gt(Mix88, return.alleles=FALSE)
colnames(geno) <- sapply(strsplit(colnames(geno), "_"), function(x) x[[1]][1])
dim(geno)
```

```
## [1] 26128714       88
```

```r
genoMix88 <- geno[as.numeric(rownames(SNP)),] 
dim(genoMix88)
```

```
## [1] 23 88
```

Save genotypes for Mix88


```r
SNP.Mix88 <- SNP
save(Mix88, genoMix88, SNP.Mix88, file="Mix88.Rdata")
write.table(genoMix88, file="genotypesMix88.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(SNP.Mix88, file="snp_data_Mix88.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
```

### Read VCF files into R: Quarter Horses


```r
QH <- read.vcfR(vcf.files["QH"], verbose = FALSE)
```

Extract VCF INFO


```r
info.snp <- data.frame(getFIX(QH))
info.snp$POS <- as.numeric(as.character(info.snp$POS))
info.snp$CHROM <- paste("chr", info.snp$CHROM, sep="")
nrow(info.snp)
```

```
## [1] 131476
```

Extract SNPs within genes of interest


```r
SNP <- lapply(1:nrow(posSNP), function(x) 
    data.frame(info.snp[as.character(info.snp$CHROM) == as.character(posSNP$chr[x]) & 
        info.snp$POS == posSNP$EquCab2[x],]))
idx <- unlist(lapply(SNP, nrow)) > 0
SNP <- data.frame(gene=posSNP$gene[idx], do.call(rbind, SNP))
SNP
```

```
##           gene CHROM      POS          ID REF ALT QUAL FILTER
## 603       BAG3  chr1 12853565  1#12853565   C   T  999   <NA>
## 3897      LDB3  chr1 83735063  1#83735063   T   C  999   <NA>
## 97389     FLNC  chr4 83738769  4#83738769   G   A  999   <NA>
## 102062    LMNA  chr5 42041844  5#42041844   A   C  999   <NA>
## 107817 PYROXD1  chr6 47640400  6#47640400   A   G  999   <NA>
## 107819 PYROXD1  chr6 47655192  6#47655192   A   G  999   <NA>
## 107820 PYROXD1  chr6 47655207  6#47655207   A   G  999   <NA>
## 107824 PYROXD1  chr6 47661977  6#47661977   G   C  999   <NA>
## 125578    PLEC  chr9 82500927  9#82500927   C   T  999   <NA>
## 125587    PLEC  chr9 82501889  9#82501889   A   G  999   <NA>
## 35208   SQSTM1 chr14  2659718  14#2659718   G   A  999   <NA>
## 37241     MYOT chr14 38519183 14#38519183   A   G  999   <NA>
```

Extract Genotypes


```r
geno <- extract.gt(QH, return.alleles=FALSE)
colnames(geno) <- sapply(strsplit(colnames(geno), "_"), function(x) x[[1]][1])
dim(geno)
```

```
## [1] 131476     39
```

```r
genoQH <- geno[as.numeric(rownames(SNP)),] 
dim(genoQH)
```

```
## [1] 12 39
```

Save genotypes for Mix88


```r
SNP.QH <- SNP
save(QH, genoQH, SNP.QH, file="QH.Rdata")
write.table(genoQH, file="genotypesQH.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(SNP.QH, file="snp_data_QH.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
```

### Run R Script


```r
# htmlRunR
# GenotypesPublicData.R nodes=1,cpus-per-task=1,time=04:00:00,mem=200G +Genotypes Public Data
```

