#!/usr/bin/env Rscript
# generate eisa reference for intron + exon quantification

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 3){
  stop("missing args: gtf fasta outpath")
}

gtf <- args[1]
fa <- args[2]
outpath <- args[3]

library(eisaR)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(GenomicFeatures)
library(stringr)

outdir <- file.path(outpath, "eisa")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

grl <- getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 40L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE
)

# grl may contain out of bounds regions (off chromosome entries)
chrom_sizes <- read.table(paste0(fa, ".fai"),
                          header = FALSE,
                          sep = "\t")
chrom_sizes <- chrom_sizes[, 1:2]
colnames(chrom_sizes) <- c("seqname", "size")
slens <- chrom_sizes$size[match(seqlevels(grl), chrom_sizes$seqname)]
seqinfo(grl) <- Seqinfo(seqnames = seqlevels(grl), seqlengths = slens)

grl <- trim(grl)

seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = Rsamtools::FaFile(fa), 
  transcripts = grl
)

# trim polyAs 
seq_names <- names(seqs)
seqs <- str_remove(as.character(seqs), "A+$")
seqs <- as(seqs, "DNAStringSet")
names(seqs) <- seq_names

mcols(seqs)$is_duplicate <- duplicated(seqs)
mcols(seqs)$seq_len <- width(seqs)
mcols(seqs)$seq_id <- match(seqs, unique(seqs))
mcols(seqs)$representative <- unlist(lapply(split(names(seqs),
                                                  mcols(seqs)$seq_id),
                                            function(x) rep(x[1], length(x))),
                                     use.names = FALSE)

seqs_no_tails_long <- seqs[width(seqs) > 25]

duplicate_seqs <- seqs_no_tails_long[mcols(seqs_no_tails_long)$is_duplicate]
duplicate_seqs <- as.data.frame(mcols(duplicate_seqs))
duplicate_seqs$DuplicateRef <- rownames(duplicate_seqs)
duplicate_seqs$RetainedRef <-  duplicate_seqs$representative
duplicate_seqs <- duplicate_seqs[, c("DuplicateRef", 
                                                "RetainedRef",
                                                "seq_len")]
rownames(duplicate_seqs) <- NULL

polished_seqs <- seqs_no_tails_long[!mcols(seqs_no_tails_long)$is_duplicate]

exportToGtf(
  grl, 
  filepath = file.path(outdir, "eisaR.gtf"))

df <- getTx2Gene(grl)

write.table(df, file.path(outdir,"tx2gene.tsv"), sep = "\t", 
            quote = FALSE, row.names = FALSE)

# write to top-level directory, also
writeXStringSet(polished_seqs, file.path(outpath,'eisa.fa'))
indexFa(file.path(outpath,'eisa.fa'))

writeXStringSet(polished_seqs, file.path(outdir,'eisa.fa'))
indexFa(file.path(outdir,'eisa.fa'))

write.table(duplicate_seqs, file.path(outdir,"duplicated_seqs.tsv"), sep = "\t", 
          quote = FALSE, row.names = FALSE)
