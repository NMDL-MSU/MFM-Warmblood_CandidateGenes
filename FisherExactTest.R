#' ### Description
#' Performe Fisher exact test on called cSNP from RNA-Seq.
#'
#' ***
#' **Code:**
#' Parent Directory:
#'
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest
#'
#' File:
#'
#' &nbsp;&nbsp;&nbsp;&nbsp;FisherExactTest.R
#'
#' **Input files:**
#' Directories
#'
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/Filtered_rst/*
#'
#' **Output files:**
#'
#' Directory:
#'
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest
#'
#' Files:
#'
#' >&nbsp;&nbsp;&nbsp;&nbsp;FisherTestResults.Rdata
#' >&nbsp;&nbsp;&nbsp;&nbsp;FisherTestResults.txt
#'
#' Render R Script
#'
#' > &nbsp;&nbsp;&nbsp;&nbsp;FisherExactTest.qsub
#'
#' ***

#' ### R Environment

#' Load required libraries
library(rcompanion)
library(viridis)
library(qvalue)

#' Clear Environment
rm(list=ls())

#' ### Load annotated cSNP
rst <- list()
dir <- "/mnt/research/NMDL/2019_WB_MFM/cSNP"
for (i in 0:5){
    load(paste(paste(dir, "/Filtered_rst", i, "/", sep=""), "Filtered_rst", i, ".Rdata", sep=""))
}

#' Merge results to list
rst <- list(rst0, rst1, rst2, rst3, rst4, rst5)

#' Genotype Information
GenoInfo <- do.call(rbind, lapply(rst, function(x) x$Info))
dim(GenoInfo)

#' Genotype Matrix
Geno <- do.call(rbind, lapply(rst, function(x) x$Geno))
dim(Geno)

#' Check that we have no missing genotypes: should be zero
sum(is.na(Geno))

#' ### Number of genotypes per chromosome
chr <- table(GenoInfo$CHROM)
names(chr) <- gsub("chr", "", names(chr))

#' Plot number of coding SNP found per chromosome.
#+ cSNP_Chr, fig.align='center', fig.width=20, fig.height=7, dpi=300
par(mar=c(5,5,4,2) + 0.1)
pl <- barplot(chr, col=(viridis(32)), ylim=c(0, 6000),
    cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5,
    xlab = "Chromosome", ylab="Coding SNP")
text(pl, chr + 100, labels=chr)

#' Prepare genotpes for fisher exact test.
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


#' ### Annotation
Ann <- lapply(rst, function(x) x$Ann)
high <- do.call(rbind, lapply(Ann, function(x)
    do.call(rbind, lapply(x, function(y) y[y$Annotation_Impact == "HIGH",]))))
nrow(high)

#' Check the number of cSNP with high impact
lk <- sort(table(high$Annotation))
lk <- lk[lk > 0]
lk

# Genes with cSNP annotated as high impact
GenesHighImpact <- unique(high$Gene_Name)
length(GenesHighImpact)

#' Table of cSNP annotation
lk <- unlist(lapply(Ann, function(x)
    unlist(lapply(x, function(y) y$Annotation))))
sort(table(lk))

#' ### Fisher Exact Test

#' Fisher Test
fisherTest <- lapply(dataM, function(x) fisher.test(x))
pval <- do.call(rbind, lapply(fisherTest, function(x) x$p.value))
colnames(pval) <- "p-value"

#' Histogram of p-values
#+ histogram_pvalues, fig.align='center', fig.width=7, fig.hight=7, dpi=600
hist(-log10(pval), col="gray")

#' Bonferoni multiple test correction: Number of significant snp
length(pval[pval[,1] < 0.1/nrow(Geno),])

#' FDR multiple test correction: Assuming the proportion of true null hypothesis equals zero (pi0=1)
qval <- qvalue(p=pval[,1])
summary(qval)

#' Screen for potential cSNP for futher evaluation.
# cSNP with p-value < 0.001
idx <- names(pval[pval < 0.001,])
idx

# Annotation of selected cSNP for screening
screenSNP <- lapply(Ann, function(x) x[names(x) %in% idx])
screenSNP <- do.call(rbind, lapply(screenSNP, function(x) do.call(rbind,x)))

#' Genes annotated to the selected SNP
info <- lapply(idx, function(x) screenSNP[grep(x, screenSNP$SNP), 1:9])
names(info) <- idx

#' ### Summary Results Fisher Exact Test
gene <- do.call(rbind, lapply(info, function(x) data.frame(
    Annotation=paste(unique(x$Annotation), collapse=","),
    Annotation_Impact=unique(x$Annotation_Impact),
    Gene_Name=paste(unique(x$Gene_Name), collapse=","),
    pvalue=fisherTest[[as.character(unique(x$SNP))]]$p.value,
    CI.L=fisherTest[[as.character(unique(x$SNP))]]$conf.int[1],
    CI.R=fisherTest[[as.character(unique(x$SNP))]]$conf.int[2],
    OddsRatio=fisherTest[[as.character(unique(x$SNP))]]$estimate)))
gene

#' Add cSNP positions
gene <- cbind(GenoInfo[idx,1:4], gene)
gene

#' Allele counts
cnt <- do.call(rbind, (lapply(idx, function(x) cbind(dataM[[x]][1,], dataM[[x]][2,]))))
rownames(cnt) <- idx
colnames(cnt) <- c("Cnt.Ref", "Cnt.Alt", "MFM.Ref", "MFM.Alt")
cnt

#' Genotype Counts
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

#' ### Save results file
# Save Rdata object
save(mfm, rst, Geno, GenoInfo, Ann, high, GenesHighImpact,
    fisherTest, qval, screenSNP, info, gene, cnt, cnt2, 
    file=paste(getwd(), "FisherTestResults.Rdata", sep="/"))

# Save results to txt file
write.table(gene, file=paste(getwd(), "FisherTestResults.txt", sep="/"),
    quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)


#' ### Run R Script

#+ eval = FALSE
# htmlRunR
# FisherExactTest.R nodes=1,cpus-per-task=1,time=02:00:00,mem=50G +Fisher Exact Test cSNPs


#' ### Session Information
sessionInfo()

