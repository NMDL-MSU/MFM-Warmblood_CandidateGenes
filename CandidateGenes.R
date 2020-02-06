
#' ### Description  
#' Extract called cSNP and Fisher Exact p-values for candidate genes.  
#' Run variant effect prediction on called cSNP for select genes.     
#'
#' ***  
#' **Code:**  
#' Parent Directory:  
#'
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest/CandidateGenes  
#'
#' File:  
#'
#' &nbsp;&nbsp;&nbsp;&nbsp;CandidateGenes.R  
#'
#' **Input files:**  
#' Directory/File  
#'
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest/FisherTestResults.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/VariantCalling/gene_location.txt  
#'
#' **Output files:**  
#'
#' Directory:  
#'
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest/CandidateGenes  
#'
#' Files:  
#'
#' >&nbsp;&nbsp;&nbsp;&nbsp;CandidateGenes.Rdata  
#'
#' Render R Script  
#'
#' > &nbsp;&nbsp;&nbsp;&nbsp;CandidateGenes.qsub  
#'
#' ***

#' Clear environment
rm(list=ls())

#' Load  Fisher Exact Test Results
dir <- "/mnt/research/NMDL/2019_WB_MFM/cSNP/FisherExactTest"
fl <- "FisherTestResults.Rdata"
load(paste(dir, fl, sep="/"))

#' Gene Positions
genePos <- read.table("/mnt/research/NMDL/2019_WB_MFM/VariantCalling/gene_location.txt")
colnames(genePos) <- c("Gene", "Chr", "ID", "Start", "End")

#' Extract cSNP for genes for interest
GenoInfo$CHROM <- gsub("chr", "", GenoInfo$CHROM)
GenoInfo$POS <- as.numeric(as.character(GenoInfo$POS))
SubSet <- lapply(1:nrow(genePos), function(x) GenoInfo[GenoInfo$CHROM %in% genePos$Chr[x] 
    & GenoInfo$POS <= genePos$End[x] & GenoInfo$POS >= genePos$Start[x],])
names(SubSet) <- genePos$Gene

#' Number of called SNP per gene
unlist(lapply(SubSet, nrow))

#' Extract Fisher Exact Test P-Values
SubSetFisher <- lapply(SubSet, function(x) unlist(lapply(fisherTest[rownames(x)], function(y) y$p.value)))
idx <- unlist(lapply(SubSetFisher, function(x) sum(x < 0.05)))
idx

#' ### Results Fisher Test Candidate Genes

#' MYOT had six cSNP with p-value < 0.05
idx <- SubSetFisher[idx > 0][[1]][SubSetFisher[idx > 0][[1]] < 0.05]
idx

#' SNP Information
GenoInfo[names(idx),c(1:6,10:23)]

#' Reduced Annotation for SNP in MYOT
do.call(rbind, lapply(Ann[[4]][names(Ann[[4]]) %in% names(idx)], function(x) x[2,1:11]))

#' Full annotation for SNPs in MYOT with Fisher Exact  p-value < 0.05
MYOT <- Ann[[4]][names(Ann[[4]]) %in% names(idx)]
names(MYOT) <- names(idx)

#' Allele counts for MYOT cSNP
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

#' Genotype counts for MYOT cSNP
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

#' MYOT cSNP Counts
CountsMYOT <- list(Allele=Cnt1, Geno=Cnt2)
CountsMYOT


#' ### Variant Effect Prediction  

#' VEP text file: SNP in candidate genes  
idx <- unlist(lapply(SubSetFisher, names))
VEP <- data.frame(chr=GenoInfo[idx, "CHROM"], 
    start=GenoInfo[idx, "POS"], end=GenoInfo[idx, "POS"], 
    allele=paste(GenoInfo[idx,"REF"], GenoInfo[idx,"ALT"], sep="/"), 
    strand=rep("+", length(idx)), identifyer=idx)
write.table(VEP, file="VEP_CandidateGenes.txt", row.names=FALSE, 
    col.names=FALSE, quote=FALSE, sep="\t")

#' Save script for VEP to file
vep <- "/mnt/home/velezdeb/R/x86_64-pc-linux-gnu-library/3.5/ensemblVEP/ensembl-vep"
fl <- "VEP_CandidateGenes.txt"
out <- "VEP_CandidateGenes_output.txt"
vep <- paste(vep, "/vep -i ", fl, " -o ", 
    out, " --cache --merged --offline --species equus_caballus", sep="")
write.table(vep, "VEP.sh", row.names=FALSE, col.names=FALSE, quote=FALSE)

#' Run VEP
#+ message = FALSE
system("bash VEP.sh")

#' Save R data
save(genePos, SubSet, SubSetFisher, MYOT, CountsMYOT, VEP, file="CandidateGenes.Rdata")

#' ### Run R Script
#+ eval = FALSE
~/bin/RunR CandidateGenes.R nodes=1,cpus-per-task=1,time=01:00:00,mem=10G +Candidate Genes cSNP


