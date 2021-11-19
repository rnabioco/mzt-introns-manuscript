
twoBitToFa \
    http://hgdownload.soe.ucsc.edu/goldenPath/xenTro9/bigZips/xenTro9.2bit \
    xenTro9.fa

samtools faidx xenTro9.fa
# 2018-07-17

rsync -avzP \
    rsync://hgdownload.cse.ucsc.edu/goldenPath/xenTro9/database/refGene.txt.gz \
    .

zcat refGene.txt.gz \
    | cut -f2- \
    | genePredToGtf -source=xenTro9.refGene.ucsc file stdin xenTro9.gtf

rm refGene.txt.gz
