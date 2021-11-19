#! /usr/bin/env bash
#BSUB -n 1
#BSUB -J dl
#BSUB -e err.txt
#BSUB -o out.txt

# download fastqs from GEO


libs=$(awk '$8 == "RNASeq" && $5 ~ /Zygote|Cleavage|Gastrula|Blastula/ {print $5":"$2}' elife-30860-supp1-v1.tsv ) 
for fq in $libs
do
    echo $fq
    sample_name=$(echo $fq | cut -d ":" -f 3) 
    tissue=$(echo $fq | cut -d ":" -f 2)
    echo $tissue
    echo $sample_name
    sampid=$(grep $sample_name  PRJEB12982.txt | cut -f 2 | uniq) 
    fqs=$(grep $sample_name  PRJEB12982.txt | cut -f 11 | uniq) 
    mkdir -p $tissue"/"$sample_name
    echo $sample_name
    
    for file_id in $fqs
    do 
      R1=$(echo $file_id | cut -d ";" -f 1)
      R2=$(echo $file_id | cut -d ";" -f 2)
      echo $R1
      echo $R2
      wget -P $tissue"/"$sample_name $R1  
      wget -P $tissue"/"$sample_name $R2
    done
done

bash format_data.sh
