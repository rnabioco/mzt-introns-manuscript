#wget http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv3.genes.gff3.gz
#gunzip Montipora_capitata_HIv3.genes.gff3.gz
#
#wget http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv3.assembly.fasta.gz
#gunzip Montipora_capitata_HIv3.assembly.fasta.gz
#
#samtools faidx Montipora_capitata_HIv3.assembly.fasta

wget http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv2.assembly.fasta.gz
gunzip Montipora_capitata_HIv2.assembly.fasta.gz
samtools faidx Montipora_capitata_HIv2.assembly.fasta

wget http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv2.genes.gff3.gz
gunzip Montipora_capitata_HIv2.genes.gff3.gz
