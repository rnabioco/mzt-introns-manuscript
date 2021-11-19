
# get genome
wget ftp://ftp.ensembl.org/pub/release-94/fasta/gallus_gallus/dna/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa.gz
gunzip Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa.gz
samtools faidx Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa

# get annotations
wget ftp://ftp.ensembl.org/pub/release-94/gtf/gallus_gallus/Gallus_gallus.Gallus_gallus-5.0.94.gtf.gz
gunzip Gallus_gallus.Gallus_gallus-5.0.94.gtf.gz
