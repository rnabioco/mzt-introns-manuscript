# ! /usr/bin/env Rscript
# gff to gtf
library(rtracklayer)
library(GenomicFeatures)

gff_fn <- "Montipora_capitata_HIv2.genes.gff3"
gff <- import(gff_fn)
txdb <- GenomicFeatures::makeTxDbFromGFF(gff_fn, format = "gff")

tx <- transcripts(txdb, use.names = FALSE)
tx$type <- "transcript"
tx$transcript_id <- tx$tx_id
tx$transcript_name <- tx$tx_name

## Identify all exons
exons <- exonsBy(txdb, use.names = TRUE)
exons <- unlist(exons)
exons$transcript_id <- as.numeric(factor(names(exons)))
exons$transcript_name <- names(exons)
exons$type <- "exon"
exons <- unname(exons)

gtf_tx <- c(tx, exons)
gtf_tx <- sort(gtf_tx)
gtf_tx$tx_id <- NULL
gtf_tx$tx_name <- NULL
gtf_tx$gene_id <- gtf_tx$transcript_id
gtf_tx$gene_name <- gtf_tx$transcript_name


export(gtf_tx, 
       "Montipora_capitata_HIv2.genes.gtf",
       format = 'gtf')
