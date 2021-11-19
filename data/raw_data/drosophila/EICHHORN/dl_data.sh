#! /usr/bin/env bash
#BSUB -n 1
#BSUB -J dl
#BSUB -e err.txt
#BSUB -o out.txt

# download fastqs from GEO


libs=$(grep "RNA-Seq" GSE83616_run_info_geo.txt | grep "wt"  | cut -f 9 )

for fq in $libs
do
    fastq-dump  --gzip $fq 
done


