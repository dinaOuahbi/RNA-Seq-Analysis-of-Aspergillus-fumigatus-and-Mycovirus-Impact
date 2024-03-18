#!/bin/bash
# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


SECONDS=0

# change working directory
cd /shared/projects/mycovirus/

# LODING ...
module load trimmomatic/0.39
module load samtools/1.15.1

echo "===> RUN TRIMMOMATIC <==="

#run trimmomatic to trim reads with poor quality 
if [ ! -d "trimmed_data/" ]; then
    mkdir "trimmed_data/"
    
else
    echo "trimmed_data/ already exists. Skipping..."
fi

# Définissez vos deux listes
liste1=($(ls reads_mycovirus/*R1*))
liste2=($(ls reads_mycovirus/*R2*))

# Vérifiez si les deux listes ont la même longueur
if [ ${#liste1[@]} -ne ${#liste2[@]} ]; then
    echo "Les listes n'ont pas la même longueur."
    exit 1
fi

# Utilisez une boucle for avec un index pour parcourir les listes simultanément
for ((i = 0; i < ${#liste1[@]}; i++)); do
    r1="${liste1[i]%%.*}"
    r2="${liste2[i]%%.*}"

    # Faites quelque chose avec les éléments des deux listes
    echo "r1 : $r1"
    echo "r2 : $r2"
    trimmomatic PE -threads 6 $r1.fastq $r2.fastq $r1+p.fastq $r1+u.fastq $r2+p.fastq $r2+u.fastq SLIDINGWINDOW:4:25
    rm $r1.fastq
    rm $r2.fastq
		
done

mkdir trimmed_data/paired/
mkdir trimmed_data/unpaired/
for file in `find reads_mycovirus/ -type f -name "*+p*"`; do echo $file;mv $file trimmed_data/paired ; done
for file in `find reads_mycovirus/ -type f -name "*+u*"`; do echo $file;mv $file trimmed_data/unpaired ; done


echo "TRIMMOMATIC DONE"
