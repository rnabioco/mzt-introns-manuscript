#!/usr/bin/env Rscript
# generate eisa masked intron reference for intron + exon quantification

info <- paste0(c("missing args: eisa_dir outpath bws ",
                 "  eisa_dir: directory containing eisaR.gtf, eisa.fa, and duplicated_seqs.tsv files",
                 "  outpath: output directory",
                 "  bws: bigwig files, supplied as space separated values.",
                 "       Note that these should be all the forward bigwigs if the" ,
                 "       library is stranded. If the library is unstranded then",
                 "       both forward and reverse bigwigs should be supplied"),
               sep = "\n")
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 3){
  stop(info)
}

eisa_dir <- args[1]
outpath <- args[2]
bws <- args[3:length(args)]

library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(GenomicFeatures)
library(Biostrings)
library(Rsamtools)

read_bigwig <- function(path, set_strand = "+") {
  bw <- rtracklayer::import(path) 
  strand(bw) <- set_strand
  bw
}

dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outpath, "intron_mask"), recursive = TRUE, showWarnings = FALSE)

# get fasta of tx seqs
fa_fn <- file.path(eisa_dir, "eisa.fa")
eisa_seqs <- import(fa_fn, format = "fasta")

# get coverage and find seqs to mask ----------------------------------------------

bw_ivls <- lapply(bws, read_bigwig)
bw_ivls <- lapply(bw_ivls, function(x) {
  seqlevels(x) <- seqlevels(eisa_seqs)
  seqinfo(x) <- seqinfo(eisa_seqs)
  x
})
bw_ivls <- unlist(as(bw_ivls, "GRangesList"))

# make disjoint intervals, summing coverage in duplicated regions
bw_cov <- disjoin(bw_ivls, with.revmap=TRUE)
mcols(bw_cov) <- aggregate(bw_ivls, mcols(bw_cov)$revmap, score = sum(score),drop = FALSE)

# select intron regions to mask
masked_ivls <- reduce(bw_cov[bw_cov$score > 0])
masked_intron_ivls <- masked_ivls[str_detect(seqnames(masked_ivls), "-I[0-9]*$")]

export(masked_intron_ivls,
       file.path(outpath, "intron_mask", "intron_mask_txcoords_eisa_bt2_align.bed"))

dup_tx <- read.table(file.path(eisa_dir, "duplicated_seqs.tsv"), header = TRUE)
gtf_fn <- file.path(eisa_dir, "eisaR.gtf")
txdb <- makeTxDbFromGFF(gtf_fn)
tx <- transcripts(txdb, use.names = TRUE)
exons <- exons(txdb, use.names = TRUE)
exons <- exons[!str_detect(names(exons), "-I[0-9]*")]
exons <- reduce(unname(exons))

masked_genomic_intron_ivls <- reduce(mapFromTranscripts(masked_intron_ivls, tx))
all_genomic_intron_ivls <- reduce(unname(tx[str_detect(names(tx), "-I[0-9]*"), ]))

export(masked_genomic_intron_ivls,
       file.path(outpath, "intron_mask", "intron_mask_genomiccoords_eisa_bt2_align.bed"))

export(tx[str_detect(names(tx), "-I[0-9]*"), ],
       file.path(outpath, "intron_mask", "eisa_introns.bed"))


# Report amount of intron masking -----------------------------------------

# - Report based on amount of genomic intronic sequence
# - Exclude exon flank (+/- 40bp) from intron seq
# - Report masking due to intron regions that overlap  exons from quant (masked due to exon overlap)

# remove exon seqs (and removes flanks)
hits <- findOverlaps(all_genomic_intron_ivls, exons)
toSubtract <- reduce(extractList(exons, as(hits, "List")), 
                     ignore.strand = TRUE)
ans <- psetdiff(all_genomic_intron_ivls, toSubtract)
all_genomic_intron_ivls <- unlist(ans[width(ans) > 0L])

ans <- findOverlapPairs(masked_genomic_intron_ivls, all_genomic_intron_ivls)
mcols(ans)$overlap_width <- width(pintersect(ans, ignore.strand = TRUE))
masked_nts <- sum(mcols(ans)$overlap_width)

prop_masked <- signif(sum(mcols(ans)$overlap_width) / sum(width(all_genomic_intron_ivls)), 3)

message("proportion of genomic intron nucleotides masked to Ns: ", prop_masked )

# Filter and mask sequences ------------------------------------------------

eisa_intron_seqs <- eisa_seqs[str_detect(names(eisa_seqs), "-I[0-9]*")]
intron_tx_coord <- GRanges(seqnames = names(eisa_intron_seqs), 
                           ranges = IRanges(start = 1, 
                                            end = seqlengths(eisa_intron_seqs)),
                           strand = "+")
seqlevels(intron_tx_coord) <- seqlevels(eisa_intron_seqs)
seqinfo(intron_tx_coord) <- seqinfo(eisa_intron_seqs) 

ans <- findOverlapPairs(masked_intron_ivls, intron_tx_coord)
mcols(ans)$overlap_width <- width(pintersect(ans, ignore.strand = TRUE))
masked_nts <- unlist(lapply(split(mcols(ans)$overlap_width, seqnames(second(ans))),
                            sum))

tx_lens <- width(intron_tx_coord)
names(tx_lens) <- seqnames(intron_tx_coord)

tx_info <- data.frame(masked_nts = masked_nts, 
           tx_len = tx_lens[names(masked_nts)],
           tx = names(masked_nts))
tx_info <- tx_info[!tx_info$tx %in% dup_tx$DuplicateRef, ]

# annotate # of Ns from genome (not due to coverage masking)
genome_N_nts <- letterFrequency(eisa_seqs, "N")[, 1]
names(genome_N_nts) <- names(eisa_seqs)

stopifnot(all(tx_info$tx %in% names(genome_N_nts)))

tx_info$genome_masked <- genome_N_nts[tx_info$tx]
tx_info$prop_masked <- (tx_info$genome_masked + tx_info$masked_nts) / tx_info$tx_len
tx_info$unmasked_nts <- tx_info$tx_len - tx_info$genome_masked - tx_info$masked_nts

# set cutoffs for introns to include for quant
to_drop <- rownames(tx_info)[tx_info$prop_masked >= 0.80 | tx_info$unmasked_nts < 25]

message("removing ", 
        length(to_drop), " (", signif(length(to_drop) / nrow(tx_info), 3), ")",
        " introns with:\n >= 80% masked sequence or < 25 nt unmasked sequence")

prop_masked <- signif(sum(tx_info$masked_nts) / sum(tx_info$tx_len), 3)

message("proportion of all intron nucleotides masked to Ns: ", prop_masked )

intron_expr <- tx_info[tx_info$masked_nts > 0, ]
prop_expr_masked <- signif(sum(intron_expr$masked_nts) / sum(intron_expr$tx_len), 6)

message("proportion of intron nucleotides masked to Ns\nfor introns with at least 1 read: ",
        prop_expr_masked )

# mask intron sequence -----------------------------------------------------------------

irl_to_mask <- as(masked_intron_ivls, "IntegerRangesList")
irl_to_mask <- irl_to_mask[elementNROWS(irl_to_mask) > 0]
seqs_to_mask <- eisa_seqs[names(irl_to_mask)]

mask_seqs <- as(strrep("N", width(masked_intron_ivls)), "DNAStringSet")
mask_seqs <- split(mask_seqs, seqnames(masked_intron_ivls))
mask_seqs <- mask_seqs[names(irl_to_mask)]

stopifnot(identical(names(mask_seqs), names(irl_to_mask)))
stopifnot(identical(names(mask_seqs), names(seqs_to_mask)))

# names not used in replaceAt
irl_to_mask <- unname(irl_to_mask)
masked_seqs <- replaceAt(seqs_to_mask, irl_to_mask, mask_seqs)

# add back masked sequences
eisa_masked_seqs <- eisa_seqs[setdiff(names(eisa_seqs), names(masked_seqs))]
eisa_masked_seqs <- c(eisa_masked_seqs, masked_seqs)

# drop duplicate and invalid introns, due to excessive masking, or being duplicate sequence
eisa_masked_seqs <- eisa_masked_seqs[!names(eisa_masked_seqs) %in% to_drop]
eisa_masked_seqs <- eisa_masked_seqs[!names(eisa_masked_seqs) %in% dup_tx$DuplicateRef]

outfa <- file.path(outpath, "eisa_masked.fa")
writeXStringSet(eisa_masked_seqs, outfa)
idx_fn <- indexFa(outfa)
