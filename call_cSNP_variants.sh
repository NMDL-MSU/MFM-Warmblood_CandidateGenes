#==============================================================================
#   File: call_cSNP_variants.sh
#   Directory code: /mnt/research/NMDL/2019_WB_MFM/cSNP
#   Date: April 22, 2019
#   Description: Run mpileup and bcftools on bam files to call all cSNPs for 
#          Warmblood MFM project 
#   Run: bash call_cSNP_variants.sh
#------------------------------------------------------------------------------
#   Input files:
#       /mnt/research/NMDL/2019_WB_MFM/Tophat/Sorted/*.bam
#       /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/DNA/Index_Chr_VariantCalling/index_vector.txt
#       /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/DNA/Index_Chr_VariantCalling/index_call_variants.txt
#
#   Output files to directory:
#       /mnt/research/NMDL/2019_WB_MFM/cSNP
#==============================================================================

### Directories
# Work directory
cd /mnt/research/NMDL/2019_WB_MFM
mkdir cSNP
dir=/mnt/research/NMDL/2019_WB_MFM/cSNP

# Move bash scripts and gene index to work directory
mv call_cSNP_variants.sh $dir

# Copy genome indexes to work directory
scp /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/DNA/Index_Chr_VariantCalling/index*.txt $dir

# Qstat directory
mkdir $dir/qstat
qstat=$dir/qstat

### Input Files
# Bam files directory
Bam=/mnt/research/NMDL/2019_WB_MFM/Tophat/Sorted

# Reference EquCab3 (same as the reference used for generating BAM files)
ref=/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq/EquCab3/bowtie2/Index_EquCab3.fa

# snpEff directory
snpEff=/mnt/research/NMDL/snpEff_latest_core/snpEff

# Chromosome indexes
chr=(`cat $dir/index_vector.txt | cut -f1 -d: | uniq`)

# Submit script to HPC
cd $dir
for ((j=1; j<=${#chr[@]} ; j++ )) do
cd $dir
idx=(`cat index_call_variants.txt | cut -f$j`)

for ((i=0; i<${#idx[@]} ; i++ )) do
echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=50G
#SBATCH -J '${idx[$i]}'
#SBATCH -o '${idx[$i]}.o%j'

# Load Modules
module purge
module load powertools
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load SAMtools/1.9
module load bcftools/1.9.64
module load matplotlib/2.1.2-Python-3.6.4
module load Java/1.8.0_192
module list

# Move to output directory
cd '$Bam'

# Run mpileup
bcftools mpileup -Ou -C50 -E -Q25 -a DV,AD,ADF,ADR,SP -f '$ref' -r '${idx[$i]}' '*.bam' \
    | bcftools call -mv -Oz \
    | bcftools filter -s LowQual -g3 -G10 \
    -e '"'"'%QUAL<30'"'"' \
    -o '$dir'/'${idx[$i]}'.vcf

# Add cSNP annotation with snpEff
cd '$snpEff'
java -Xmx4g -jar snpEff.jar GCF_002863925.1_EquCab3.0 '$dir'/'${idx[$i]}'.vcf > '$dir'/'${idx[$i]}'_annot.vcf
rm '$dir'/'${idx[$i]}'.vcf

# Zip VCF file
gzip '$dir'/'${idx[$i]}'_annot.vcf

# Job details
echo 'Job Details'
scontrol show job $SLURM_JOB_ID' > $qstat/${idx[$i]}.qsub

# Submit script to hpcc
cd $qstat
sbatch ${idx[$i]}.qsub

# Move to index directory
cd $index

done
done

