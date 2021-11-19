#!/usr/bin/env Rscript
# generate eisa masked intron reference for intron + exon quantification

args = commandArgs(trailingOnly=TRUE)
 
if(length(args) < 5){
   stop("missing args: eisa_dir species outpath threads fwd_bws")
}
 
eisa_dir <- args[1]
species <- args[2]
outpath <- args[3]
n_cores <- as.integer(args[4])
fwd_bws <- args[5:length(args)]

library(purrr)
library(valr)
library(dplyr)
library(stringr)
library(rtracklayer)
library(doParallel)
library(readr)

read_bigwig <- function(path, set_strand = "+") {
      # note that rtracklayer will produce a one-based GRanges object
      rtracklayer::import(path) %>%
        dplyr::as_tibble(.) %>%
        dplyr::mutate(chrom = as.character(seqnames),
                                        start = start - 1L,
                                                          strand =
                                        set_strand) %>%
        dplyr::select(chrom, start, end, score, strand)
}

tidy_gtf <- function(gtf_fn, zero_based_coords = TRUE){
    gtf <- rtracklayer::import(gtf_fn)
    gtf <- as.data.frame(gtf)
    gtf <- dplyr::mutate_if(gtf, is.factor, as.character)
    res <- dplyr::rename(gtf, chrom = seqnames)

      if(zero_based_coords) {
              res <- dplyr::mutate(res, start = start - 1L)
      }

      res <- tibble::as_tibble(res)
        res
}

dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outpath, "intron_mask"), recursive = TRUE, showWarnings = FALSE)

fwd_df <- map_dfr(fwd_bws,
                  ~read_bigwig(.x, set_strand = "+")) 

chunk_ivls <- function(df, n) {
  if (n == 1){
    return(list(df))
  }
  chroms <- unique(df$chrom)
  chrom_splits <- split(chroms, cut(seq_along(chroms), n, labels = FALSE))
  map(chrom_splits, 
      ~filter(df, chrom %in% .x))
}

cl <- makeCluster(n_cores)  
registerDoParallel(cl)  
message("making cluster with ", n_cores, " cores")

fwd_df <- foreach(i = chunk_ivls(fwd_df, n_cores),
                  .packages = c("valr"),
                  .combine = rbind
) %dopar% {
  bed_partition(i, score = sum(score))
}

stopCluster(cl)

ivl_df <- mutate(fwd_df, strand = "+") %>% 
  dplyr::select(chrom:score, strand)

masked_ivls <- filter(ivl_df, score > 0) %>% 
  group_by(strand) %>% 
  bed_merge() %>% 
  ungroup()

gtf_fn <- file.path(eisa_dir, "eisaR.gtf")
gtf_tbl <- tidy_gtf(gtf_fn)

exon_bed <- gtf_tbl %>% 
  filter(type == "exon",
         str_detect(ID, "-I[0-9]*$")) %>% 
  dplyr::select(chrom, start, end, transcript_id, gene_id, score, strand)


mask_df <- masked_ivls %>% 
  filter(strand == "+") %>% 
  group_by(strand)

# convert to tx coords
pre_mrna_bed <- exon_bed %>%
  mutate(end = end - start,
         start = 0,
         chrom = transcript_id,
         strand = "+") %>% 
  group_by(strand)

tx_coords <- bed_intersect(mask_df, pre_mrna_bed, 
                           suffix = c("", ".y")) %>% 
  mutate(name = ".",
         score = 0) %>% 
  select(chrom, start, end, name, score, strand) %>% 
  unique()

# set some basic cutoffs for introns to include for quant
to_drop <- bed_intersect(mask_df , pre_mrna_bed, 
                         suffix = c("", ".y")) %>% 
  group_by(chrom) %>% 
  summarize(masked_nts = sum(.overlap), 
            total_nts = unique(end.y - start.y),
            prop_masked = masked_nts / total_nts, 
            unmasked_nts = total_nts - masked_nts) %>% 
  filter(prop_masked > 0.80 |  unmasked_nts < 25)

message("removing ", 
        length(unique(to_drop$chrom)),
        " introns with:\n >= 80% masked sequence\n or < 25 nt unmasked sequence")

total_intron_nts <- filter(exon_bed, strand == "+") %>%
  bed_merge() %>%
  mutate(width = end - start) %>% 
  summarize(nt_length = sum(width)) %>% 
  pull(nt_length)

masked_intron_nts <- tx_coords %>% 
  bed_merge() %>% 
  mutate(width = end - start) %>% 
  summarize(nt_length = sum(width)) %>% 
  pull(nt_length)

prop_masked <- signif(masked_intron_nts / total_intron_nts, 3)
message("proportion of intron nucleotides masked to Ns: ", prop_masked )
  
db_out <- tx_coords %>% 
  dplyr::select(chrom,
                start, 
                end, 
                name,
                score,
                strand) %>% 
  bed_sort() %>% 
  mutate_if(is.numeric, format, scientific = F)  %>% 
  mutate_all(str_trim)

write_tsv(db_out,
          file.path(outpath, "intron_mask", "intron_mask_txcoords_eisa_bt2_align.bed"),
          col_names = F)

dup_tx <- read_tsv(file.path(eisa_dir, "duplicated_seqs.tsv"))

chroms <- read_tsv(file.path(eisa_dir, "eisa.fa.fai"),
                   col_names = c("chrom", "size"),
                   col_type = c("ci---"))

chroms <- chroms %>% 
  filter(!chrom %in% dup_tx$DuplicateRef,
         !chrom %in% to_drop$chrom)

write_lines(chroms$chrom, file.path(outpath, "intron_mask","eisa_unique_transcripts.tsv"))

