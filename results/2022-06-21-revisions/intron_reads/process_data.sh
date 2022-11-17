#BSUB -n 12
#BSUB -R "select[mem>30] rusage[mem=30] span[hosts=1]"
#BSUB -J introns
#BSUB -e logs/err_%J.txt
#BSUB -o logs/out_%J.txt

mkdir -p logs
mkdir -p rissland
projdir=~/Projects/rissland-mzt-introns/mzt-introns-manuscript
bash \
  $projdir/results/2022-06-21-revisions/src/realign-intronic-reads.sh \
  $projdir/dbases/rissland/eisa.fa.fai \
  $projdir/data/rissland/bt2/eisa/ \
  $projdir/results/2022-06-21-revisions/intron_reads/rissland \
  rissland/intron_reads_ \
  $projdir/dbases/drosophila/star/BDGP6/ \
  $projdir/dbases/drosophila/ext/Drosophila_melanogaster.BDGP6.84.gtf

mkdir -p eichhorn
bash \
  $projdir/results/2022-06-21-revisions/src/realign-intronic-reads.sh \
  $projdir/dbases/eichhorn/eisa.fa.fai \
  $projdir/data/eichhorn/bt2/eisa/ \
  $projdir/results/2022-06-21-revisions/intron_reads/eichhorn \
  eichhorn/intron_reads_ \
  $projdir/dbases/drosophila/star/BDGP6/ \
  $projdir/dbases/drosophila/ext/Drosophila_melanogaster.BDGP6.84.gtf

mkdir -p white
bash \
  $projdir/results/2022-06-21-revisions/src/realign-intronic-reads-pe.sh \
  $projdir/dbases/white/eisa.fa.fai \
  $projdir/data/white/bt2/eisa/ \
  $projdir/results/2022-06-21-revisions/intron_reads/white \
  white/intron_reads_ \
  $projdir/dbases/zebrafish/star/danRer10 \
  $projdir/dbases/zebrafish/ext/danRer10.gtf

mkdir -p owens_polya
bash \
  $projdir/results/2022-06-21-revisions/src/realign-intronic-reads-pe.sh \
  $projdir/dbases/owens_polya/eisa.fa.fai \
  $projdir/data/owens_polya/bt2/eisa/ \
  $projdir/results/2022-06-21-revisions/intron_reads/owens_polya \
  owens_polya/intron_reads_ \
  $projdir/dbases/xenopus_t/star/xenTro9 \
  $projdir/dbases/xenopus_t/ext/xenTro9.gtf

mkdir -p owens_ribo
bash \
  $projdir/results/2022-06-21-revisions/src/realign-intronic-reads-pe.sh \
  $projdir/dbases/owens_ribo/eisa.fa.fai \
  $projdir/data/owens_ribo/bt2/eisa/ \
  $projdir/results/2022-06-21-revisions/intron_reads/owens_ribo \
  owens_ribo/intron_reads_ \
  $projdir/dbases/xenopus_t/star/xenTro9 \
  $projdir/dbases/xenopus_t/ext/xenTro9.gtf

mkdir -p chille
bash \
  $projdir/results/2022-06-21-revisions/src/realign-intronic-reads-pe.sh \
  $projdir/dbases/chille/eisa.fa.fai \
  $projdir/data/chille/bt2/eisa/ \
  $projdir/results/2022-06-21-revisions/intron_reads/chille \
  chille/intron_reads_ \
  $projdir/dbases/coral/star/Montipora_capitata_HIv2 \
  $projdir/dbases/coral/ext/Montipora_capitata_HIv2.genes.gtf


