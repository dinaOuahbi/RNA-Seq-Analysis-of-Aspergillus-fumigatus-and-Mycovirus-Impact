#!/bin/bash
# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


SECONDS=0

# change working directory
cd /shared/projects/mycovirus/

# LODING ...
module load fastqc/0.11.9
module load multiqc/1.13

echo "===> RUN FASTQC <==="

fastqc_res="fastqc_res/"
multiqc_res="multiqc_res/"
fastq_folder="reads_mycovirus/"

###

if [ ! -d "$fastqc_res" ]; then
    mkdir "$fastqc_res"
    
else
    echo "$fastqc_res already exists. Skipping..."
fi

###

if [ ! -d "$multiqc_res" ]; then
    mkdir "$multiqc_res"
    
else
    echo "$multiqc_res already exists. Skipping..."
fi

###

for file in `ls $fastq_folder`; do
    echo $file
    fastqc $fastq_folder$file -o $fastqc_res
done

echo "===> RUN MULTIQC <==="

multiqc -f -o $multiqc_res $fastqc_res


