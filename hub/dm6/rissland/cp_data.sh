datadir="/beevol/home/riemondy/Projects/mzt-introns/data/bigwigs/drosophila/RISSLAND"

#cp $datadir/*fwd.bw .
#cp $datadir/*rev_neg.bw .

# convert ensembl to ucsc

for bw in *.bw
do echo $bw
    bigWigToBedGraph \
        $bw stdout \
        | python3 rename_chrom.py \
          -bg - \
          -c ~/src/ChromosomeMappings/BDGP6_ensembl2UCSC.txt \
          | sort -k1,1 -k2,2n \
          > $bw".tmp"
    
    bedGraphToBigWig \
      $bw".tmp" \
      http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes \
      $bw

    rm $bw".tmp"

done
