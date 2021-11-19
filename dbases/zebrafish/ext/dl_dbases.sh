rsync -avzP \
    rsync://hgdownload.cse.ucsc.edu/goldenPath/danRer10/database/refGene.txt.gz \
    . 

zcat refGene.txt.gz \
    | cut -f2- \
    | genePredToGtf -source=danRer10.refGene.ucsc file stdin danRer10.gtf

rm refGene.txt.gz

twoBitToFa http://hgdownload.soe.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.2bit danRer10.fa
samtools faidx danRer10.fa
