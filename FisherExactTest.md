---
title: Fisher Exact Test cSNPs
author: Deborah Velez-Irizarry
date: Tue Aug 20 12:34:13 EDT 2019
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description
Performe Fisher exact test on called cSNP from RNA-Seq.

***
**Code:**
Parent Directory:

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest

File:

&nbsp;&nbsp;&nbsp;&nbsp;FisherExactTest.R

**Input files:**
Directories

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/Filtered_rst/*

**Output files:**

Directory:

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest

Files:

>&nbsp;&nbsp;&nbsp;&nbsp;FisherTestResults.Rdata
>&nbsp;&nbsp;&nbsp;&nbsp;FisherTestResults.txt

Render R Script

> &nbsp;&nbsp;&nbsp;&nbsp;FisherExactTest.qsub

***
### R Environment
Load required libraries


```r
library(rcompanion)
library(viridis)
```

```
## Loading required package: viridisLite
```

```r
library(qvalue)
```

Clear Environment


```r
rm(list=ls())
```

### Load annotated cSNP


```r
rst <- list()
dir <- "/mnt/research/NMDL/2019_WB_MFM/cSNP"
for (i in 0:5){
    load(paste(paste(dir, "/Filtered_rst", i, "/", sep=""), "Filtered_rst", i, ".Rdata", sep=""))
}
```

Merge results to list


```r
rst <- list(rst0, rst1, rst2, rst3, rst4, rst5)
```

Genotype Information


```r
GenoInfo <- do.call(rbind, lapply(rst, function(x) x$Info))
dim(GenoInfo)
```

```
## [1] 72365    27
```

Genotype Matrix


```r
Geno <- do.call(rbind, lapply(rst, function(x) x$Geno))
dim(Geno)
```

```
## [1] 72365    16
```

Check that we have no missing genotypes: should be zero


```r
sum(is.na(Geno))
```

```
## [1] 0
```

### Number of genotypes per chromosome


```r
chr <- table(GenoInfo$CHROM)
names(chr) <- gsub("chr", "", names(chr))
```

Plot number of coding SNP found per chromosome.


```r
par(mar=c(5,5,4,2) + 0.1)
pl <- barplot(chr, col=(viridis(32)), ylim=c(0, 6000),
    cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5,
    xlab = "Chromosome", ylab="Coding SNP")
text(pl, chr + 100, labels=chr)
```

<img src="figure/cSNP_Chr-1.svg" title="plot of chunk cSNP_Chr" alt="plot of chunk cSNP_Chr" style="display: block; margin: auto;" />

Prepare genotpes for fisher exact test.


```r
# Sum the number of allele for reference and alternate allele for each group (MFM, Control)
mfm <- c("9056", "9066", "12711", "12130", "12429", "9054", "12124", "12065")
dataM <- lapply(1:nrow(Geno), function(x)
    data.frame(Ref=rbind(
        Cnt=sum(Geno[x,!colnames(Geno) %in% mfm] %in% "0/0") * 2 +
        sum(Geno[x,!colnames(Geno) %in% mfm] %in% "0/1"),

        MFM=sum(Geno[x,colnames(Geno) %in% mfm] %in% "0/0") * 2 +
        sum(Geno[x,colnames(Geno) %in% mfm] %in% "0/1")),

        Alt=rbind(
        Cnt=sum(Geno[x,!colnames(Geno) %in% mfm] %in% "1/1") * 2 +
        sum(Geno[x,!colnames(Geno) %in% mfm] %in% "0/1"),

        MFM=sum(Geno[x,colnames(Geno) %in% mfm] %in% "1/1") * 2 +
        sum(Geno[x,colnames(Geno) %in% mfm] %in% "0/1"))))
names(dataM) <- rownames(Geno)
```

### Annotation


```r
Ann <- lapply(rst, function(x) x$Ann)
high <- do.call(rbind, lapply(Ann, function(x)
    do.call(rbind, lapply(x, function(y) y[y$Annotation_Impact == "HIGH",]))))
nrow(high)
```

```
## [1] 1247
```

Check the number of cSNP with high impact


```r
lk <- sort(table(high$Annotation))
lk <- lk[lk > 0]
lk
```

```
## 
##       start_lost&splice_region_variant 
##                                      2 
## splice_acceptor_variant&intron_variant 
##                                      9 
##                              stop_lost 
##                                     16 
##    splice_donor_variant&intron_variant 
##                                     18 
##                             start_lost 
##                                     48 
##      stop_gained&splice_region_variant 
##                                     84 
##                            stop_gained 
##                                   1070
```

```r
# Genes with cSNP annotated as high impact
GenesHighImpact <- unique(high$Gene_Name)
length(GenesHighImpact)
```

```
## [1] 277
```

Table of cSNP annotation


```r
lk <- unlist(lapply(Ann, function(x)
    unlist(lapply(x, function(y) y$Annotation))))
sort(table(lk))
```

```
## lk
##                         start_lost&splice_region_variant 
##                                                        2 
##                   splice_acceptor_variant&intron_variant 
##                                                        9 
##                                    stop_retained_variant 
##                                                       14 
##                                                stop_lost 
##                                                       16 
##                                  initiator_codon_variant 
##                                                       17 
##                      splice_donor_variant&intron_variant 
##                                                       18 
##                                               start_lost 
##                                                       48 
##                        stop_gained&splice_region_variant 
##                                                       84 
##                     splice_region_variant&intron_variant 
##                                                      234 
##                                    splice_region_variant 
##                                                      276 
## splice_region_variant&non_coding_transcript_exon_variant 
##                                                      280 
##                                              stop_gained 
##                                                     1070 
##           5_prime_UTR_premature_start_codon_gain_variant 
##                                                     1419 
##                   missense_variant&splice_region_variant 
##                                                     2113 
##                 splice_region_variant&synonymous_variant 
##                                                     2271 
##                                        intergenic_region 
##                                                     5424 
##                                      5_prime_UTR_variant 
##                                                     7712 
##                       non_coding_transcript_exon_variant 
##                                                    10304 
##                                    upstream_gene_variant 
##                                                    28092 
##                                         missense_variant 
##                                                    50182 
##                                           intron_variant 
##                                                    54896 
##                                  downstream_gene_variant 
##                                                    64401 
##                                       synonymous_variant 
##                                                    76372 
##                                      3_prime_UTR_variant 
##                                                    89525
```

### Fisher Exact Test
Fisher Test


```r
fisherTest <- lapply(dataM, function(x) fisher.test(x))
pval <- do.call(rbind, lapply(fisherTest, function(x) x$p.value))
colnames(pval) <- "p-value"
```

Histogram of p-values


```r
hist(-log10(pval), col="gray")
```

<img src="figure/histogram_pvalues-1.svg" title="plot of chunk histogram_pvalues" alt="plot of chunk histogram_pvalues" style="display: block; margin: auto;" />

Bonferoni multiple test correction: Number of significant snp


```r
length(pval[pval[,1] < 0.1/nrow(Geno),])
```

```
## [1] 0
```

FDR multiple test correction: Assuming the proportion of true null hypothesis equals zero (pi0=1)


```r
qval <- qvalue(p=pval[,1])
summary(qval)
```

```
## 
## Call:
## qvalue(p = pval[, 1])
## 
## pi0:	1	
## 
## Cumulative number of significant calls:
## 
##           <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
## p-value        1      8   123    359   780 1338 28959
## q-value        0      0     0      0     0    0     0
## local FDR      0      0     0      0     0    0     0
```

Screen for potential cSNP for futher evaluation.


```r
# cSNP with p-value < 0.001
idx <- names(pval[pval < 0.001,])
idx
```

```
## [1] "chr1_173870607_84316"  "chr1_173871701_84325"  "chr8_41341153_438206" 
## [4] "chr9_7230262_459580"   "chr11_5963963_553576"  "chr15_32523098_728217"
## [7] "chr17_67642251_824751" "chr20_30790745_908837"
```

```r
# Annotation of selected cSNP for screening
screenSNP <- lapply(Ann, function(x) x[names(x) %in% idx])
screenSNP <- do.call(rbind, lapply(screenSNP, function(x) do.call(rbind,x)))
```

Genes annotated to the selected SNP


```r
info <- lapply(idx, function(x) screenSNP[grep(x, screenSNP$SNP), 1:9])
names(info) <- idx
```

### Summary Results Fisher Exact Test


```r
gene <- do.call(rbind, lapply(info, function(x) data.frame(
    Annotation=paste(unique(x$Annotation), collapse=","),
    Annotation_Impact=unique(x$Annotation_Impact),
    Gene_Name=paste(unique(x$Gene_Name), collapse=","),
    pvalue=fisherTest[[as.character(unique(x$SNP))]]$p.value,
    CI.L=fisherTest[[as.character(unique(x$SNP))]]$conf.int[1],
    CI.R=fisherTest[[as.character(unique(x$SNP))]]$conf.int[2],
    OddsRatio=fisherTest[[as.character(unique(x$SNP))]]$estimate)))
gene
```

```
##                                                                              Annotation
## chr1_173870607_84316                        3_prime_UTR_variant,downstream_gene_variant
## chr1_173871701_84325                        3_prime_UTR_variant,downstream_gene_variant
## chr8_41341153_438206                          5_prime_UTR_variant,upstream_gene_variant
## chr9_7230262_459580                                                      intron_variant
## chr11_5963963_553576  3_prime_UTR_variant,upstream_gene_variant,downstream_gene_variant
## chr15_32523098_728217                                                    intron_variant
## chr17_67642251_824751        3_prime_UTR_variant,downstream_gene_variant,intron_variant
## chr20_30790745_908837                         3_prime_UTR_variant,upstream_gene_variant
##                       Annotation_Impact                 Gene_Name
## chr1_173870607_84316           MODIFIER          FAM177A1,PPP2R3C
## chr1_173871701_84325           MODIFIER          FAM177A1,PPP2R3C
## chr8_41341153_438206           MODIFIER LOC100058979,LOC102150007
## chr9_7230262_459580            MODIFIER                    ZNF704
## chr11_5963963_553576           MODIFIER RNF157,LOC111775585,FOXJ1
## chr15_32523098_728217          MODIFIER                    PAIP2B
## chr17_67642251_824751          MODIFIER               STK24,FARP1
## chr20_30790745_908837          MODIFIER               DDR1,GTF2H4
##                             pvalue     CI.L         CI.R OddsRatio
## chr1_173870607_84316  8.157212e-04 0.000000    0.3142941   0.00000
## chr1_173871701_84325  8.157212e-04 0.000000    0.3142941   0.00000
## chr8_41341153_438206  5.130768e-04 0.000000    0.3206498   0.00000
## chr9_7230262_459580   8.157212e-04 3.181733          Inf       Inf
## chr11_5963963_553576  6.351663e-04 2.974107 1495.7194678  28.71752
## chr15_32523098_728217 1.538373e-05 0.000000    0.1567439   0.00000
## chr17_67642251_824751 6.351663e-04 2.974107 1495.7194678  28.71752
## chr20_30790745_908837 9.650523e-04 2.620301  237.4962739  18.42318
```

Add cSNP positions


```r
gene <- cbind(GenoInfo[idx,1:4], gene)
gene
```

```
##                       CHROM       POS REF ALT
## chr1_173870607_84316   chr1 173870607   G   A
## chr1_173871701_84325   chr1 173871701   T   G
## chr8_41341153_438206   chr8  41341153   G   C
## chr9_7230262_459580    chr9   7230262   A   G
## chr11_5963963_553576  chr11   5963963   C   T
## chr15_32523098_728217 chr15  32523098   C   T
## chr17_67642251_824751 chr17  67642251   A   G
## chr20_30790745_908837 chr20  30790745   A   G
##                                                                              Annotation
## chr1_173870607_84316                        3_prime_UTR_variant,downstream_gene_variant
## chr1_173871701_84325                        3_prime_UTR_variant,downstream_gene_variant
## chr8_41341153_438206                          5_prime_UTR_variant,upstream_gene_variant
## chr9_7230262_459580                                                      intron_variant
## chr11_5963963_553576  3_prime_UTR_variant,upstream_gene_variant,downstream_gene_variant
## chr15_32523098_728217                                                    intron_variant
## chr17_67642251_824751        3_prime_UTR_variant,downstream_gene_variant,intron_variant
## chr20_30790745_908837                         3_prime_UTR_variant,upstream_gene_variant
##                       Annotation_Impact                 Gene_Name
## chr1_173870607_84316           MODIFIER          FAM177A1,PPP2R3C
## chr1_173871701_84325           MODIFIER          FAM177A1,PPP2R3C
## chr8_41341153_438206           MODIFIER LOC100058979,LOC102150007
## chr9_7230262_459580            MODIFIER                    ZNF704
## chr11_5963963_553576           MODIFIER RNF157,LOC111775585,FOXJ1
## chr15_32523098_728217          MODIFIER                    PAIP2B
## chr17_67642251_824751          MODIFIER               STK24,FARP1
## chr20_30790745_908837          MODIFIER               DDR1,GTF2H4
##                             pvalue     CI.L         CI.R OddsRatio
## chr1_173870607_84316  8.157212e-04 0.000000    0.3142941   0.00000
## chr1_173871701_84325  8.157212e-04 0.000000    0.3142941   0.00000
## chr8_41341153_438206  5.130768e-04 0.000000    0.3206498   0.00000
## chr9_7230262_459580   8.157212e-04 3.181733          Inf       Inf
## chr11_5963963_553576  6.351663e-04 2.974107 1495.7194678  28.71752
## chr15_32523098_728217 1.538373e-05 0.000000    0.1567439   0.00000
## chr17_67642251_824751 6.351663e-04 2.974107 1495.7194678  28.71752
## chr20_30790745_908837 9.650523e-04 2.620301  237.4962739  18.42318
```

Allele counts


```r
cnt <- do.call(rbind, (lapply(idx, function(x) cbind(dataM[[x]][1,], dataM[[x]][2,]))))
rownames(cnt) <- idx
colnames(cnt) <- c("Cnt.Ref", "Cnt.Alt", "MFM.Ref", "MFM.Alt")
cnt
```

```
##                       Cnt.Ref Cnt.Alt MFM.Ref MFM.Alt
## chr1_173870607_84316        7       9      16       0
## chr1_173871701_84325        7       9      16       0
## chr8_41341153_438206        6       8      16       0
## chr9_7230262_459580        16       0       7       9
## chr11_5963963_553576       15       1       5      11
## chr15_32523098_728217       0      12      13       3
## chr17_67642251_824751      11       5       1      15
## chr20_30790745_908837      14       2       4      12
```

Genotype Counts


```r
control <- colnames(Geno)[!colnames(Geno) %in% mfm]
cnt2 <- do.call(rbind, apply(Geno[idx,],1, function(x) data.frame(
    A=sum(x[control] == "0/0"),
    B=sum(x[control] == "0/1"),
    C=sum(x[control] == "1/1"),
    D=sum(x[mfm] == "0/0"),
    E=sum(x[mfm] == "0/1"),
    F=sum(x[mfm] == "1/1"))))
colnames(cnt2) <- c("Cnt.0/0", "Cnt.0/1", "Cnt.1/1",
    "MFM.0/0", "MFM.0/1", "MFM.1/1")
cnt2
```

```
##                       Cnt.0/0 Cnt.0/1 Cnt.1/1 MFM.0/0 MFM.0/1 MFM.1/1
## chr1_173870607_84316        1       5       2       8       0       0
## chr1_173871701_84325        1       5       2       8       0       0
## chr8_41341153_438206        2       2       3       8       0       0
## chr9_7230262_459580         8       0       0       1       5       2
## chr11_5963963_553576        7       1       0       0       5       3
## chr15_32523098_728217       0       0       6       5       3       0
## chr17_67642251_824751       4       3       1       0       1       7
## chr20_30790745_908837       6       2       0       0       4       4
```

### Save results file


```r
# Save Rdata object
save(mfm, rst, Geno, GenoInfo, Ann, high, GenesHighImpact,
    fisherTest, qval, screenSNP, info, gene, cnt, cnt2, 
    file=paste(getwd(), "FisherTestResults.Rdata", sep="/"))

# Save results to txt file
write.table(gene, file=paste(getwd(), "FisherTestResults.txt", sep="/"),
    quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
```

### Run R Script


```r
# htmlRunR
# FisherExactTest.R nodes=1,cpus-per-task=1,time=02:00:00,mem=50G +Fisher Exact Test cSNPs
```

### Session Information


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
## [1] en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] qvalue_2.14.1     viridis_0.5.1     viridisLite_0.3.0 rcompanion_2.2.1 
## [5] knitr_1.20       
## 
## loaded via a namespace (and not attached):
##  [1] zoo_1.8-2          modeltools_0.2-22  nortest_1.0-4     
##  [4] tidyselect_0.2.5   coin_1.2-2         purrr_0.3.2       
##  [7] reshape2_1.4.3     splines_3.5.1      lattice_0.20-38   
## [10] colorspace_1.3-2   expm_0.999-2       stats4_3.5.1      
## [13] survival_2.42-3    rlang_0.3.4        pillar_1.4.1      
## [16] foreign_0.8-71     glue_1.3.1         multcomp_1.4-8    
## [19] plyr_1.8.4         multcompView_0.1-7 stringr_1.4.0     
## [22] munsell_0.5.0      gtable_0.3.0       mvtnorm_1.0-10    
## [25] codetools_0.2-15   evaluate_0.12      lmtest_0.9-36     
## [28] manipulate_1.0.1   highr_0.7          TH.data_1.0-8     
## [31] Rcpp_1.0.1         scales_1.0.0       gridExtra_2.3     
## [34] ggplot2_3.2.0      stringi_1.4.3      dplyr_0.8.1       
## [37] grid_3.5.1         tools_3.5.1        sandwich_2.4-0    
## [40] magrittr_1.5       DescTools_0.99.28  lazyeval_0.2.2    
## [43] tibble_2.1.1       crayon_1.3.4       pkgconfig_2.0.2   
## [46] MASS_7.3-51.1      Matrix_1.2-14      assertthat_0.2.1  
## [49] R6_2.4.0           boot_1.3-20        EMT_1.1           
## [52] compiler_3.5.1
```

