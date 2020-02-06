### Description:  
#' Extract genotypes of interest from public data sources.  
#' [Warmblood PRJEB23301](https://www.ebi.ac.uk/ena/data/view/PRJEB23301)
#' [Mix Breeds PRJEB28306](https://www.ebi.ac.uk/ena/data/view/PRJEB28306)
#' [Quarter Horses PRJEB30116](https://www.ebi.ac.uk/ena/data/view/PRJEB30116)
#'  
#' ***  
#' **Code:**  
#' Parent Directory:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/Public_Data  
#'  
#' File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;GenotypesPublicData.R  
#'  
#' **Input files:**  
#' Directories  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/Public_Data  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/DNA/VariantCalling  
#'   
#' Files:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;379_RAO_HD_Warmbloods.ne1000.vcf  
#' >&nbsp;&nbsp;&nbsp;&nbsp;horses.88.vars.flt.pass.ebi.vcf  
#' >&nbsp;&nbsp;&nbsp;&nbsp;cSNP.vcf  
#'  
#' **Output files:**  
#'  
#' Directory:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/2019_WB_MFM/Public_Data/GenotypesPublicData  
#'  
#' Files:  
#'  
#' >&nbsp;&nbsp;&nbsp;&nbsp;TB.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;genoTB.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;dataTB.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;snp_data_TB.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;WB379.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;genoWB379.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;snp_data_WB379.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Mix88.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;genoMix88.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;snp_data_Mix88.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;QH.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;genoQH.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;snp_data_QH.txt  
#'  
#' Render R Script  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;GenotypesPublicData.qsub  
#'  
#' ***  

#' ### R Environment

#' Load required libraries
library(vcfR)

#' Clear Environment
rm(list=ls())

#' Session Information
sessionInfo()


#' ### Load Data

#' VCF cSNP File
Dir1 <- "/mnt/research/NMDL/2019_WB_MFM/Public_Data"
Dir2 <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/DNA/VariantCalling"
vcf.files <- c(paste(Dir1, list.files(Dir1)[grep("vcf", list.files(Dir1))], sep="/"),
    paste(Dir2, list.files(Dir2)[grep("cSNP.vcf", list.files(Dir2))], sep="/"))
names(vcf.files) <- c("WB379", "QH", "Mix88", "TB")
vcf.files


#' SNP of interest
posSNP <- read.table("/mnt/research/NMDL/2019_WB_MFM/Public_Data/snp_of_interest.txt", header=TRUE)

#' ### Read VCF files into R: Thoroughbred Horses
TB <- read.vcfR(vcf.files["TB"], verbose = FALSE)

#' Extract VCF INFO
info.snp <- data.frame(getFIX(TB))
info.snp$POS <- as.numeric(as.character(info.snp$POS))
nrow(info.snp)

#' Extract SNPs within genes of interest
SNP <- lapply(1:nrow(posSNP), function(x) 
    data.frame(gene=posSNP$gene[x], info.snp[as.character(info.snp$CHROM) == as.character(posSNP$chr[x]) & 
        info.snp$POS == posSNP$EquCab3[x],]))
SNP <- do.call(rbind, SNP)
SNP

#' Extract Genotypes
geno <- extract.gt(TB, return.alleles=FALSE)
colnames(geno) <- sapply(strsplit(colnames(geno), "_"), function(x) x[[1]][1])
dim(geno)
genoTB <- geno[as.numeric(rownames(SNP)),] 
dim(genoTB)

#' Animal Information RER Thoroughbred
RER <- read.table("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq/RER_Thoroughbred_Animal_Information.txt", 
    header=TRUE, sep="\t", row.names=1)
RER$ID <- as.character(RER$ID)
rownames(RER) <- RER$ID

#' Animal Information for Glycogen Thoroughbred
Gly <- read.table("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/Glycogen_Project_Information.txt",
    header=TRUE, sep="\t")

# Retain information on sequenced animals 
Gly <- Gly[!is.na(Gly$MSMS_Plate),]
Gly$TimePoint <- factor(as.character(Gly$TimePoint), exclude="Rep48h", 
    levels=c("Pre", "Depl", "Rep24h", "Rep72h"))
Gly$Diet <- factor(as.character(Gly$Diet), exclude="DietGoldenMax")
rownames(Gly) <- paste("G", Gly$MSMS_ID, sep="")
head(Gly)

# Animal Ids for Glycogen Thoroughbreds
idx <- lapply(unique(Gly$Animal), function(x) rownames(Gly[Gly$Animal == x,]))
names(idx) <- unique(Gly$Animal)
idx

#' Check that genotypes are consistent across replicates
ck <- lapply(idx, function(x) unlist(apply(genoTB[,x], 1, function(y) sum(!y %in% y[1]))))
lapply(ck, function(x) x[x>0])

#' Remove repeated animals from genotype matrix
un <- unlist(lapply(idx, function(x) x[1]))
unG <-  genoTB[,(ncol(geno)-40):ncol(geno)][,un]
colnames(unG) <- names(un)

#' RER genotypes
rerG <- genoTB[,colnames(genoTB) %in% RER$ID]

#' TB Genotype Matrix
genoTB <- cbind(unG, rerG)
dim(genoTB)

#' SNP Info for TB
SNPinfo.TB <- SNP


#' Data matrix for TB
dataTB <- rbind(data.frame(ID=Gly[un, "Animal"], Dx=rep("Control", length(un)), Breed=rep("Thoroughbred", length(un))),
    data.frame(RER[colnames(rerG),c("ID", "Dx")], Breed=rep("Thoroughbred", ncol(rerG))))
rownames(dataTB) <- NULL

#' Save genotypes for Thoroughbreds
SNP.TB <- SNP
save(TB, genoTB, SNP.TB, dataTB, file="TB.Rdata")
write.table(genoTB, file="genotypesTB.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(SNP.TB, file="snp_data_TB.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(dataTB, file="dataTB.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

#' ### Read VCF files into R: Warmblood Horses 379
WB379 <- read.vcfR(vcf.files["WB379"], verbose = FALSE)

#' Extract VCF INFO
info.snp <- data.frame(getFIX(WB379))
info.snp$POS <- as.numeric(as.character(info.snp$POS))
nrow(info.snp)

#' Extract SNPs within genes of interest
SNP <- lapply(1:nrow(posSNP), function(x) 
    data.frame(info.snp[as.character(info.snp$CHROM) == as.character(posSNP$chr[x]) & 
        info.snp$POS == posSNP$EquCab2[x],]))
idx <- unlist(lapply(SNP, nrow)) > 0
SNP <- data.frame(gene=posSNP$gene[idx], do.call(rbind, SNP))
SNP

#' Extract Genotypes
geno <- extract.gt(WB379, return.alleles=FALSE)
colnames(geno) <- sapply(strsplit(colnames(geno), "_"), function(x) x[[1]][1])
dim(geno)
genoWB379 <- geno[as.numeric(rownames(SNP)),] 
dim(genoWB379)

#' Save genotypes for Warmblood
SNP.WB379 <- SNP
save(WB379, genoWB379, SNP.WB379, file="WB379.Rdata")
write.table(genoWB379, file="genotypesWB379.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(SNP.WB379, file="snp_data_WB379.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")


#' ### Read VCF files into R: Mixed Breed Horses 88
system.time(
    Mix88 <- read.vcfR(vcf.files["Mix88"], verbose = FALSE)
)

#' Extract VCF INFO
info.snp <- data.frame(getFIX(Mix88))
info.snp$POS <- as.numeric(as.character(info.snp$POS))
nrow(info.snp)

#' Extract SNPs within genes of interest
SNP <- lapply(1:nrow(posSNP), function(x) 
    data.frame(info.snp[as.character(info.snp$CHROM) == as.character(posSNP$chr[x]) & 
        info.snp$POS == posSNP$EquCab3[x],]))
idx <- unlist(lapply(SNP, nrow)) > 0
SNP <- data.frame(gene=posSNP$gene[idx], do.call(rbind, SNP))
SNP

#' Extract Genotypes
geno <- extract.gt(Mix88, return.alleles=FALSE)
colnames(geno) <- sapply(strsplit(colnames(geno), "_"), function(x) x[[1]][1])
dim(geno)
genoMix88 <- geno[as.numeric(rownames(SNP)),] 
dim(genoMix88)

#' Save genotypes for Mix88
SNP.Mix88 <- SNP
save(Mix88, genoMix88, SNP.Mix88, file="Mix88.Rdata")
write.table(genoMix88, file="genotypesMix88.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(SNP.Mix88, file="snp_data_Mix88.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")


#' ### Read VCF files into R: Quarter Horses
QH <- read.vcfR(vcf.files["QH"], verbose = FALSE)

#' Extract VCF INFO
info.snp <- data.frame(getFIX(QH))
info.snp$POS <- as.numeric(as.character(info.snp$POS))
info.snp$CHROM <- paste("chr", info.snp$CHROM, sep="")
nrow(info.snp)

#' Extract SNPs within genes of interest
SNP <- lapply(1:nrow(posSNP), function(x) 
    data.frame(info.snp[as.character(info.snp$CHROM) == as.character(posSNP$chr[x]) & 
        info.snp$POS == posSNP$EquCab2[x],]))
idx <- unlist(lapply(SNP, nrow)) > 0
SNP <- data.frame(gene=posSNP$gene[idx], do.call(rbind, SNP))
SNP

#' Extract Genotypes
geno <- extract.gt(QH, return.alleles=FALSE)
colnames(geno) <- sapply(strsplit(colnames(geno), "_"), function(x) x[[1]][1])
dim(geno)
genoQH <- geno[as.numeric(rownames(SNP)),] 
dim(genoQH)

#' Save genotypes for Mix88
SNP.QH <- SNP
save(QH, genoQH, SNP.QH, file="QH.Rdata")
write.table(genoQH, file="genotypesQH.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(SNP.QH, file="snp_data_QH.txt", 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

#' ### Run R Script
#+ eval = FALSE
# htmlRunR
# GenotypesPublicData.R nodes=1,cpus-per-task=1,time=04:00:00,mem=200G +Genotypes Public Data

