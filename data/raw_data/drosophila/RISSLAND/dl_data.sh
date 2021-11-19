#! /usr/bin/env bash

# download fastqs from GEO
libs=$(awk 'BEGIN {FS="\t"} $10 == "none" && $11 == "Embryo" && $12 != "png50" {print $6}' GSE98106_run_info_geo.txt )

for fq in $libs
do
    fastq-dump  --gzip $fq 
done

