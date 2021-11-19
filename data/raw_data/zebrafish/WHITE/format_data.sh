#! /usr/bin/env bash
#BSUB -n 1
#BSUB -J mv 
#BSUB -e err.txt
#BSUB -o out.txt


for d in */*/;
do echo $d; 
    echo "combining fastqs"

    fqs1=$(echo $d*_1.fastq.gz)
    fqs2=$(echo $d*_2.fastq.gz)

    outdir=$(dirname $d)
    prefix=$(basename $d)

    cat $fqs1 >  $prefix"_1.fastq.gz"
    cat $fqs2 >  $prefix"_2.fastq.gz"
done
