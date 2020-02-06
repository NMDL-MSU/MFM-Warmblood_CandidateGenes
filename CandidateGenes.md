---
title: Candidate Genes cSNP
author: Deborah Velez-Irizarry
date: Mon Oct 7 18:19:42 EDT 2019
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description  
Extract called cSNP and Fisher Exact p-values for candidate genes.  
Run variant effect prediction on called cSNP for select genes.     

***  
**Code:**  
Parent Directory:  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest/CandidateGenes  

File:  

&nbsp;&nbsp;&nbsp;&nbsp;CandidateGenes.R  

**Input files:**  
Directory/File  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest/FisherTestResults.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/VariantCalling/gene_location.txt  

**Output files:**  

Directory:  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest/CandidateGenes  

Files:  

>&nbsp;&nbsp;&nbsp;&nbsp;CandidateGenes.Rdata  

Render R Script  

> &nbsp;&nbsp;&nbsp;&nbsp;CandidateGenes.qsub  

***
Clear environment


```r
rm(list=ls())
```

Load  Fisher Exact Test Results


```r
dir <- "/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest"
fl <- "FisherTestResults.Rdata"
load(paste(dir, fl, sep="/"))
```

Gene Positions


```r
genePos <- read.table("/mnt/research/NMDL/2019_WB_MFM/VariantCalling/gene_location.txt")
colnames(genePos) <- c("Gene", "Chr", "ID", "Start", "End")
```

Extract cSNP for genes for interest


```r
GenoInfo$CHROM <- gsub("chr", "", GenoInfo$CHROM)
GenoInfo$POS <- as.numeric(as.character(GenoInfo$POS))
SubSet <- lapply(1:nrow(genePos), function(x) GenoInfo[GenoInfo$CHROM %in% genePos$Chr[x] 
    & GenoInfo$POS <= genePos$End[x] & GenoInfo$POS >= genePos$Start[x],])
names(SubSet) <- genePos$Gene
```

Number of called SNP per gene


```r
unlist(lapply(SubSet, nrow))
```

```
##     Des   CRYAB    MYOT    ZASP    FLNC    BAG3    FHL1  DNAJB6    PLEC 
##       8       1      58     126      15      10       5      52      40 
##    LMNA   ACTA1   HSPB8      KY PYROXD1   SQSTM    TIA1 
##       6      12      33      36       9      10       5
```

Extract Fisher Exact Test P-Values


```r
SubSetFisher <- lapply(SubSet, function(x) unlist(lapply(fisherTest[rownames(x)], function(y) y$p.value)))
idx <- unlist(lapply(SubSetFisher, function(x) sum(x < 0.05)))
idx
```

```
##     Des   CRYAB    MYOT    ZASP    FLNC    BAG3    FHL1  DNAJB6    PLEC 
##       0       0       6       0       0       0       0       0       0 
##    LMNA   ACTA1   HSPB8      KY PYROXD1   SQSTM    TIA1 
##       0       0       0       0       0       0       0
```

### Results Fisher Test Candidate Genes
MYOT had six cSNP with p-value < 0.05


```r
idx <- SubSetFisher[idx > 0][[1]][SubSetFisher[idx > 0][[1]] < 0.05]
idx
```

```
## chr14_37819145_690635 chr14_37822768_690653 chr14_37823712_690654 
##            0.04338154            0.04338154            0.04338154 
## chr14_37824086_690659 chr14_37827558_690673 chr14_37827935_690679 
##            0.04338154            0.04338154            0.04338154
```

SNP Information


```r
GenoInfo[names(idx),c(1:6,10:23)]
```

```
##                       CHROM      POS REF ALT  QUAL FILTER   DP       VDB
## chr14_37819145_690635    14 37819145   T   C 999.0   PASS 1110  0.171437
## chr14_37822768_690653    14 37822768   G   A 999.0   PASS  594  0.240978
## chr14_37823712_690654    14 37823712   G   A 999.0   PASS  541  0.271348
## chr14_37824086_690659    14 37824086   T   C 999.0   PASS  483 0.0799765
## chr14_37827558_690673    14 37827558   A   T 795.0   PASS  307  0.468762
## chr14_37827935_690679    14 37827935   A   T 821.0   PASS  287 0.0743581
##                            RPB         MQB      BQB     MQSB     SGB MQ0F
## chr14_37819145_690635 0.284807 1.44222e-32 0.990584 0.974519 61.4493    0
## chr14_37822768_690653 0.425153 1.22039e-34  0.61885 0.826989  25.524    0
## chr14_37823712_690654 0.553241 1.01163e-21 0.723567 0.976965 23.2062    0
## chr14_37824086_690659 0.139554 7.81694e-33 0.982547 0.953288 55.4547    0
## chr14_37827558_690673 0.692527 8.64578e-17 0.924013 0.904021 8.72015    0
## chr14_37827935_690679 0.749218 1.37203e-17  0.96431 0.380087 27.2508    0
##                            ICB       HOB AC AN           DP4 MQ
## chr14_37819145_690635 0.149204 0.0488281  5 32 482,368,77,50 48
## chr14_37822768_690653 0.149204 0.0488281  5 32 262,182,30,35 49
## chr14_37823712_690654 0.149204 0.0488281  5 32 238,152,37,24 48
## chr14_37824086_690659 0.149204 0.0488281  5 32 204,131,41,25 49
## chr14_37827558_690673 0.149204 0.0488281  5 32   133,91,28,8 49
## chr14_37827935_690679 0.149204 0.0488281  5 32   133,61,30,8 48
```

Reduced Annotation for SNP in MYOT


```r
do.call(rbind, lapply(Ann[[4]][names(Ann[[4]]) %in% names(idx)], function(x) x[2,1:11]))
```

```
##                                         SNP Allele     Annotation
## chr14_37819145_690635 chr14_37819145_690635      C intron_variant
## chr14_37822768_690653 chr14_37822768_690653      A intron_variant
## chr14_37823712_690654 chr14_37823712_690654      A intron_variant
## chr14_37824086_690659 chr14_37824086_690659      C intron_variant
## chr14_37827558_690673 chr14_37827558_690673      T intron_variant
## chr14_37827935_690679 chr14_37827935_690679      T intron_variant
##                       Annotation_Impact Gene_Name   Gene_ID Feature_Type
## chr14_37819145_690635          MODIFIER      MYOT gene19424   transcript
## chr14_37822768_690653          MODIFIER      MYOT gene19424   transcript
## chr14_37823712_690654          MODIFIER      MYOT gene19424   transcript
## chr14_37824086_690659          MODIFIER      MYOT gene19424   transcript
## chr14_37827558_690673          MODIFIER      MYOT gene19424   transcript
## chr14_37827935_690679          MODIFIER      MYOT gene19424   transcript
##                           Feature_ID Transcript_BioType Rank
## chr14_37819145_690635 XM_023617646.1     protein_coding  5/9
## chr14_37822768_690653 XM_014730661.2     protein_coding  2/9
## chr14_37823712_690654 XM_014730661.2     protein_coding  2/9
## chr14_37824086_690659 XM_014730661.2     protein_coding  2/9
## chr14_37827558_690673 XM_023617646.1     protein_coding  1/9
## chr14_37827935_690679 XM_023617646.1     protein_coding  1/9
##                               HGVS.c
## chr14_37819145_690635   c.131+271A>G
## chr14_37822768_690653   c.357-367C>T
## chr14_37823712_690654  c.357-1311C>T
## chr14_37824086_690659  c.357-1685A>G
## chr14_37827558_690673  c.-281-647T>A
## chr14_37827935_690679 c.-281-1024T>A
```

Full annotation for SNPs in MYOT with Fisher Exact  p-value < 0.05


```r
MYOT <- Ann[[4]][names(Ann[[4]]) %in% names(idx)]
names(MYOT) <- names(idx)
```

Allele counts for MYOT cSNP


```r
cont <- colnames(Geno)[!colnames(Geno) %in% mfm]
Cnt1 <- list()
for(i in names(idx)){
    x <- Geno[i,cont]
    one <- data.frame(Cont.Ref=sum(x == "0/0")*2 + sum(x == "0/1"),
        Cont.Alt=sum(x == "1/1")*2 + sum(x == "0/1"))

    x <- Geno[i,mfm]
    two <- data.frame(MFM.Ref=sum(x == "0/0")*2 + sum(x == "0/1"),
        MFM.Alt=sum(x == "1/1")*2 + sum(x == "0/1"))

    Cnt1[[i]] <- data.frame(one, two)
}
Cnt1 <- do.call(rbind, Cnt1)
```

Genotype counts for MYOT cSNP


```r
cont <- colnames(Geno)[!colnames(Geno) %in% mfm]
Cnt2 <- list()
for(i in names(idx)){
    x <- Geno[i,cont]
    one <- data.frame(Cont.Ref=sum(x == "0/0"), Cont.Het=sum(x == "0/1"),
        Cont.Alt=sum(x == "1/1"))

    x <- Geno[i,mfm]
    two <- data.frame(MFM.Ref=sum(x == "0/0"), MFM.Het= sum(x == "0/1"),
        MFM.Alt=sum(x == "1/1"))

    Cnt2[[i]] <- data.frame(one, two)
}
Cnt2 <- do.call(rbind, Cnt2)
```

MYOT cSNP Counts


```r
CountsMYOT <- list(Allele=Cnt1, Geno=Cnt2)
CountsMYOT
```

```
## $Allele
##                       Cont.Ref Cont.Alt MFM.Ref MFM.Alt
## chr14_37819145_690635       11        5      16       0
## chr14_37822768_690653       11        5      16       0
## chr14_37823712_690654       11        5      16       0
## chr14_37824086_690659       11        5      16       0
## chr14_37827558_690673       11        5      16       0
## chr14_37827935_690679       11        5      16       0
## 
## $Geno
##                       Cont.Ref Cont.Het Cont.Alt MFM.Ref MFM.Het MFM.Alt
## chr14_37819145_690635        3        5        0       8       0       0
## chr14_37822768_690653        3        5        0       8       0       0
## chr14_37823712_690654        3        5        0       8       0       0
## chr14_37824086_690659        3        5        0       8       0       0
## chr14_37827558_690673        3        5        0       8       0       0
## chr14_37827935_690679        3        5        0       8       0       0
```

### Variant Effect Prediction  
VEP text file: SNP in candidate genes  


```r
idx <- unlist(lapply(SubSetFisher, names))
VEP <- data.frame(chr=GenoInfo[idx, "CHROM"], 
    start=GenoInfo[idx, "POS"], end=GenoInfo[idx, "POS"], 
    allele=paste(GenoInfo[idx,"REF"], GenoInfo[idx,"ALT"], sep="/"), 
    strand=rep("+", length(idx)), identifyer=idx)
write.table(VEP, file="VEP_CandidateGenes.txt", row.names=FALSE, 
    col.names=FALSE, quote=FALSE, sep="\t")
```

Save script for VEP to file


```r
vep <- "/mnt/home/velezdeb/R/x86_64-pc-linux-gnu-library/3.5/ensemblVEP/ensembl-vep"
fl <- "VEP_CandidateGenes.txt"
out <- "VEP_CandidateGenes_output.txt"
vep <- paste(vep, "/vep -i ", fl, " -o ", 
    out, " --cache --merged --offline --species equus_caballus", sep="")
write.table(vep, "VEP.sh", row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Run VEP


```r
system("bash VEP.sh")
```

Save R data


```r
save(genePos, SubSet, SubSetFisher, MYOT, CountsMYOT, VEP, file="CandidateGenes.Rdata")
```

### Run R Script


```r
~/bin/RunR CandidateGenes.R nodes=1,cpus-per-task=1,time=01:00:00,mem=10G +Candidate Genes cSNP
```

