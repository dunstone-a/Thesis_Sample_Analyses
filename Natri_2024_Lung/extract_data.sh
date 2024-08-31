#!/bin/bash
# Extract Natri GEO data 
# Amelia Dunstone
# 2024-05-15

cd ~/Scratch/Natri_2024_Lung/Millie_Analysis/data

cd raw

# Extract tar files in place
for f in *tar
do 
  echo ${f}
  tar -xf ${f}
done


cd ../limix

# Extract tar.gz files in place
for f in *tar.gz
do 
  echo ${f}
  tar -xzf ${f}
done

ls *.txt > filenames_limix
head endothelial_aCap_qtl_results_all.txt 

cd ../mashr

# Extract mashr files
for f in *tar.gz
do 
  echo ${f}
  tar -xzf ${f}
done

ls *.tsv > filenames_mashr
head endothelial_aCap.tsv

# Extract .gz seurat file to .rds
cd ../seurat

for f in *.gz
do 
    echo ${f}
    gunzip ${f}
done




 
