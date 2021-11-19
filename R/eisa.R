#!/usr/bin/env Rscript
# generate eisa reference for intron + exon quantification

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 4){
  stop("missing args: gtf fasta species outpath")
}

gtf <- args[1]
fa <- args[2]
species <- args[3]
outpath <- args[4]

library(eisaR)
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)

outdir <- file.path(outpath, "eisa")
dir.create(outdir, recursive = TRUE, showWarnings = TRUE)

grl <- getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 40L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE
)

# grl may contain out of bounds regions (off chromosome entries)
chrom_sizes <- read_tsv(str_c(fa, ".fai"),
         col_types = 'ci---', 
         col_names = c("chrom", "size"))

chrom_sizes <- left_join(data.frame(chrom = seqlevels(grl)), chrom_sizes)
seqinfo(grl) <- Seqinfo(seqnames = chrom_sizes$chrom, seqlengths = chrom_sizes$size)

grl <- trim(grl)

seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = Rsamtools::FaFile(fa), 
  transcripts = grl
)
# trim polyAs 
seqs_no_tails <- map(as.character(seqs), ~str_remove(.x, "A+$")) %>% unlist()

seqs_no_tails <- tibble(
  id = names(seqs_no_tails),
  seq = seqs_no_tails
)

seqs_no_tails <- seqs_no_tails %>%
  group_by(seq) %>%
  mutate(representative = dplyr::first(id), 
         is_duplicate = representative != id) %>% 
  ungroup() %>% 
  mutate(seq_len = nchar(seq))

too_short <- filter(seqs_no_tails, seq_len < 25)
seqs_no_tails_long <- filter(seqs_no_tails, seq_len > 25)

duplicate_seqs <- filter(seqs_no_tails_long, is_duplicate) %>% 
  dplyr::select(RetainedRef = representative,
                DuplicateRef = id,
                everything())

polished_seqs <- filter(seqs_no_tails_long, !is_duplicate) 

seq_strings <- polished_seqs$seq
names(seq_strings) <- polished_seqs$id

seq_strings <- DNAStringSet(seq_strings)

exportToGtf(
  grl, 
  filepath = file.path(outdir, "eisaR.gtf"))

df <- getTx2Gene(grl)

write_tsv(df, file.path(outdir,"tx2gene.tsv"))

# write to top-level directory and eisa directory
writeXStringSet(seq_strings, file.path(outpath,'eisa.fa'))
indexFa(file.path(outpath,'eisa.fa'))

writeXStringSet(seq_strings, file.path(outdir,'eisa.fa'))
indexFa(file.path(outdir,'eisa.fa'))

write_tsv(duplicate_seqs, file.path(outdir,"duplicated_seqs.tsv"))
