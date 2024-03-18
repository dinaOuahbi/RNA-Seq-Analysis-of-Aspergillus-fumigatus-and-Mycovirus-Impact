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
module load trimmomatic/0.39
module load samtools/1.15.1
module load hisat2/2.2.1
module load htseq/0.13.5

echo ""
echo "RUN FASTQC"
echo ""

#mkdir fastqc_res/
#for file in `ls data/`
#do
#fastqc data/$file -o fastqc_res/
#done

echo ""
echo "RUN MULTIQC"
echo ""

#mkdir multiqc_res/
#multiqc -f -o multiqc_res/ fastqc_res/

echo ""
echo "RUN TRIMMOMATIC"
echo ""

#run trimmomatic to trim reads with poor quality 

#mkdir trimmed_data
# Définissez vos deux listes
#liste1=($(ls data/*R1*))
#liste2=($(ls data/*R2*))

# Vérifiez si les deux listes ont la même longueur
#if [ ${#liste1[@]} -ne ${#liste2[@]} ]; then
    #echo "Les listes n'ont pas la même longueur."
#    exit 1
#fi

# Utilisez une boucle for avec un index pour parcourir les listes simultanément
#for ((i = 0; i < ${#liste1[@]}; i++)); do
    #r1="${liste1[i]%%.*}"
    #r2="${liste2[i]%%.*}"

    # Faites quelque chose avec les éléments des deux listes
    #echo "r1 : $r1"
    #echo "r2 : $r2"

#    #trimmomatic PE -threads 6 $r1.fastq $r2.fastq $r1+p.fastq $r1+u.fastq $r2+p.fastq $r2+u.fastq SLIDINGWINDOW:4:25

#done
#mkdir trimmed_data
#for file in `find data/ -type f -name "*+p*"`; do echo $file;mv $file trimmed_data/paired ; done
#for file in `find data/ -type f -name "*+u*"`; do echo $file;mv $file trimmed_data/unpaired ; done

#echo "TRIMMOMATIC DONE"


#trimmomatic SE -threads 4 data/demo.fastq data/demo_trimmed.fastq TRAILING:10 -phred33
#echo "Trimmomatic finished running!"


#fastqc data/demo_trimmed.fastq -o data/


echo ""
echo "RUN HISAT2"
ECHO ""


# run alignment
#mkdir HISAT2/data
# Définissez vos deux listes
#liste1=($(ls trimmed_data/paired/*R1*))
#liste2=($(ls trimmed_data/paired/*R2*))

# Vérifiez si les deux listes ont la même longueur
#if [ ${#liste1[@]} -ne ${#liste2[@]} ]; then
    #echo "Les listes n'ont pas la même longueur."
    #exit 1
#fi

# Utilisez une boucle for avec un index pour parcourir les listes simultanément
#for ((i = 0; i < ${#liste1[@]}; i++)); do
    #r1="${liste1[i]%%.*}"
    #r2="${liste2[i]%%.*}"

    # Faites quelque chose avec les éléments des deux listes
    #echo "r1 : $r1"
    #echo "r2 : $r2"

    #r1_base=$(basename $r1)
    #r2_base=$(basename $r2)

    #hisat2 -p 12 --dta -x HISAT2/af293/af293 -1 $r1.fastq -2 $r2.fastq | samtools sort -o HISAT2/data/$r1_base.bam
#done



#hisat2 -q --rna-strandness R -x HISAT2/af293/af293 -U data/demo_trimmed.fastq | samtools sort -o HISAT2/demo_trimmed.bam
#echo "HISAT2 finished running!" 

#samtools view -h demo_trimmed.bam | less

echo ""
echo "RUN HTSEQ-COUNT"
echo ""

#wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-57/gff3/aspergillus_fumigatus/Aspergillus_fumigatus.ASM265v1.57.gff3.gz
#gunzip Aspergillus_fumigatus.ASM265v1.57.gff3.gz


#chr1    GenBank gene    1000    9000    .   +   .   ID=gene001;Name=GeneA

# quantification

for bamfile in `ls HISAT2/data`; do
    echo "====> $bamfile runing <===="
    htseq-count -f bam HISAT2/data/$bamfile HISAT2/Aspergillus_fumigatus.ASM265v1.57.gtf > results/counts.txt
done
echo "Quantification with htseq-count is finished"





