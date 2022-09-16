
#BSUB -n 12
#BSUB -R "select[mem>30] rusage[mem=30] span[hosts=1]"
#BSUB -J introns
#BSUB -e logs/err_%J.txt
#BSUB -o logs/out_%J.txt


eisa_fai=$1
bamdir=$2
outdir=$3
outid=$4
star_idx=$5
gtf=$6
sif=~/Projects/rissland-mzt-introns/mzt-introns-manuscript/docker/mzt-introns.sif

sample=$(basename $2 | sed 's/.bam//')
intron_pos=$(dirname $1)"/intron_tx_coords.bed"
awk 'BEGIN{OFS=FS="\t"} $1 ~ /-I/ {print $1, 40, $2 - 40}' $1 \
   | awk 'BEGIN{OFS=FS="\t"} $3 > $2' > $intron_pos

for b in $bamdir/*.bam
do echo $b
  sample=$(basename $b | sed 's/.bam//')
  singularity exec $sif \
    samtools view -b -L $intron_pos $b \
    | samtools fastq \
    | gzip > $outdir"/"$sample".fastq.gz" 
done

fqs=$outdir"/*.fastq.gz"
fqstr=$(echo $fqs | tr " " ",")

singularity exec $sif \
  STAR \
  --genomeDir $star_idx  \
  --runThreadN 12 \
  --readFilesIn $fqstr \
  --sjdbGTFfile $gtf \
  --readFilesCommand gunzip -c \
  --outFileNamePrefix $outid \
  --outFilterMultimapNmax 200 \
  --winAnchorMultimapNmax 200 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --outSAMmultNmax 200 \
  --outMultimapperOrder Random \
  --outSAMprimaryFlag AllBestScore \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --alignSJoverhangMin 8  \
  --alignSJDBoverhangMin 1 \
  --chimSegmentMin 20 \
  --chimOutType WithinBAM \
  --peOverlapNbasesMin 10 \
  --peOverlapMMp 0.1 \
  --sjdbScore 1 \
  --outFilterType Normal \
  --outSAMtype BAM Unsorted \
  --outSAMmode Full \
  --limitSjdbInsertNsj=1500000 \
  --outSAMattributes All  \
  --outSAMattrIHstart 0  \
  --outSAMstrandField intronMotif

