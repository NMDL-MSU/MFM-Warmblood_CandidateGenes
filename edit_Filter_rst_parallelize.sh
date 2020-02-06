#===================================================================
#   File: edit_Filter_rst_parallelize.sh
#   Directory code: /mnt/research/NMDL/2019_WB_MFM/cSNP  
#   Date: May 17, 2019
#   Description: Edit Filtered_rst.R to run in parallel. 
#   Run: bash edit_Filter_rst_parallelize.sh
#-------------------------------------------------------------------
#   Input directory/file:
#       /mnt/research/NMDL/2019_WB_MFM/cSNP/Filtered_rst.R
#
#   Output file to directory:
#       /mnt/research/NMDL/2019_WB_MFM/cSNP
#
#   Output file:
#      Index_cSNP.txt
#      Filtered_rst0.R
#      Filtered_rst1.R
#      Filtered_rst2.R
#===================================================================

## Directories
# Work directory
wk=/mnt/research/NMDL/2019_WB_MFM/cSNP

# Filter R script
FiltR=$wk/Filtered_rst.R

# Indexes
Rscript $wk/Index_cSNP.sh
Start=(`cat $wk/Index_cSNP.txt | cut -f1`)
End=(`cat $wk/Index_cSNP.txt | cut -f2`)

## Edit Filter_rst.R Script
cd $wk
for ((i=0; i<${#Start[@]} ; i++ )) do
    sed 's/index/'${Start[$i]}':'${End[$i]}'/g' Filtered_rst.R > $i.R
    sed 's/rst/rst'$i'/g' $i.R > Filtered_rst$i.R
    rm $i.R
done

