
#==============================================================================
#   File: check_call_variants.sh
#   Directory code: /mnt/research/NMDL/2019_WB_MFM/cSNP
#   Date: 04/22/2019
#   Description: Check mpileup bcftools and snpEff finished with output and then merge
#       vcf files
#   Run: bash check_call_variants.sh
#------------------------------------------------------------------------------
#   Input files in directory:
#       /mnt/research/NMDL/2019_WB_MFM/cSNP
#
#   Output files to directory:
#       /mnt/research/NMDL/2019_WB_MFM/cSNP
#
#   Output:
#       index_chr.txt
#       cSNP.vcf.gz
#==============================================================================

# Work Directory
out=/mnt/research/NMDL/2019_WB_MFM/cSNP

# Move script to work directory
scp -p /mnt/research/NMDL/2019_WB_MFM/check_call_variants.sh $out

# Check that you have output for each index
cd $out
idx=(`cat index_vector.txt`)

for ((i=0; i<${#idx[@]} ; i++ )) do
ls $out/${idx[$i]}_annot.vcf.gz > temp
done
rm temp

# Concatenate by chromosome
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load VCFtools/0.1.15-Perl-5.26.1

# Chromosomes 1:31 & X
cat $out/index_vector.txt > $out/index_chr.txt
sed -i 's/$/_annot.vcf.gz/' $out/index_chr.txt 

cd $out
vcf-concat -f $out/index_chr.txt  | gzip -c > cSNP.vcf.gz
mv cSNP.vcf.gz $out
rm Rm
