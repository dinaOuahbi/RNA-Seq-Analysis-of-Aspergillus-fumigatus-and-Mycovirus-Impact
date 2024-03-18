#!/bin/bash
# python 3
# author : https://github.com/dinaOuahbi
# date : 12032024
# project : mycovirus


SECONDS=0

# change working directory
cd /shared/projects/mycovirus/

# LODING ...
module load samtools/1.15.1
module load hisat2/2.2.1
module load htseq/0.13.5
module load subread/2.0.1


#fastqc data/demo_trimmed.fastq -o data/
echo "===> RUN HISAT2 <==="

# CREATE DIR
if [ ! -d "HISAT2/data/" ]; then
    mkdir "HISAT2/data/"
    
else
    echo "HISAT2/data/ already exists. Skipping..."
fi

# DEFINE TWO LISTS
liste1=($(ls trimmed_data/paired/*R1*))
liste2=($(ls trimmed_data/paired/*R2*))

# VERIFY IF LISTS HAVE SAME LENGTH
if [ ${#liste1[@]} -ne ${#liste2[@]} ]; then
    echo "Les listes n'ont pas la même longueur."
    exit 1
fi

# RUN ALIGNEMENT ON EACH SOUCHE (r1+r2)
for ((i = 0; i < ${#liste1[@]}; i++)); do
    r1="${liste1[i]%%.*}"
    r2="${liste2[i]%%.*}"

    # Faites quelque chose avec les éléments des deux listes
    echo "r1 : $r1"
    echo "r2 : $r2"

    r1_base=$(basename $r1)
    r2_base=$(basename $r2)

    hisat2 -p 12 --dta -x af_indexs/af -1 $r1.fastq -2 $r2.fastq | samtools sort -o HISAT2/data/$r1_base.bam
done



#samtools view -h demo_trimmed.bam | less

echo "===> RUN HTSEQ COUNT <==="

#wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-57/gff3/aspergillus_fumigatus/Aspergillus_fumigatus.ASM265v1.57.gff3.gz
#gunzip Aspergillus_fumigatus.ASM265v1.57.gff3.gz

#chr1    GenBank gene    1000    9000    .   +   .   ID=gene001;Name=GeneA

for bamfile in `ls HISAT2/data`; do
    echo "====> $bamfile runing <===="
    filename=$(basename $bamfile)
    featureCounts -S 2 -a Aspergillus_fumigatus.ASM265v1.57.gtf -o results/$filename.txt HISAT2/data/$bamfile
done
echo "Quantification with htseq-count is finished"





