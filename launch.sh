#!/bin/bash

#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 6
#SBATCH --mem 5GB


Rscript script/2-clustering/s1-kmean_clustering.R 
Rscript script/2-clustering/s2-agglo_clustering.R

Rscript script/3-DA/s1-DA.R/
python script/3-DA/s2-rename_genes.py

python script/4-PEA/s1bis_go_names.py
Rscript script/4-PEA/s1-gaf_to_gmt.R
Rscript script/4-PEA/s2-gsea_res.R
Rscript script/4-PEA/s3-gsea_plots.R
python script/4-PEA/s4-stack_pdf.py


echo "######################### DONE ! #######################################"
