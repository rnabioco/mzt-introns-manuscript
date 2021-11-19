#! /usr/bin/env bash

eisa=$1
bed=$2
duptx=$3

bedtools maskfasta \
    -fi $eisa \
    -bed $bed \
    -fo ${eisa/.fa/_masked.fa.tmp}

# exclude non-masked duplicate transcripts
samtools faidx ${eisa/.fa/_masked.fa.tmp} 

samtools faidx \
    -r $duptx \
    ${eisa/.fa/_masked.fa.tmp} \
    > ${eisa/.fa/_masked.fa}

samtools faidx ${eisa/.fa/_masked.fa} 

rm ${eisa/.fa/_masked.fa.tmp} 
rm ${eisa/.fa/_masked.fa.tmp.fai}

