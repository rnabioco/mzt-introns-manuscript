# globals shared across markdown docs

library(readr)
library(valr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(Matrix)
library(purrr)
library(viridis)
library(jsonlite)
library(rtracklayer)
library(ComplexHeatmap)
library(tximport) 
library(DESeq2)
library(tibble)
library(here)
library(ggrepel)
library(readxl)
library(eisaR)
library(doParallel)


#### Paths ####

project_dir <- here()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
docs_dir <- file.path(project_dir, "docs")
db_dir <- file.path(project_dir, "dbases")

#### functions ####
#' Read in a bigwig file into a valr compatible bed tibble
#' @description This function will output a 5 column tibble with
#' zero-based chrom, start, end, score, and strand columns.
#' 
#' @param path path to bigWig file
#' @param set_strand strand to add to output (defaults to "+")
#' @export 
read_bigwig <- function(path, set_strand = "+") {
  # note that rtracklayer will produce a one-based GRanges object
  rtracklayer::import(path) %>% 
    dplyr::as_tibble(.) %>% 
    dplyr::mutate(chrom = as.character(seqnames),
                  start = start - 1L, 
                  strand = set_strand) %>% 
    dplyr::select(chrom, start, end, score, strand)
}

#' Convert GTF from rtracklayer into tidy bed format
#' @param gtf_fn path to gtf file
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

#' When writing out excel workbooks using openxlsx::write.xlsx()
#' this function will set the class attributes for a column, which
#' enforces a column type in the resulting xlsx file. 
#' Useful for avoid gene names being clobbered to dates and 
#' setting scientific number formatting

set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}

#' Add transcript and gene annotations to gtf
#' 
#' 
#' 
add_tx_annotations <- function(gtf_df){
  
  if(!all(c("transcript_id", "gene_id") %in% colnames(gtf))){
    stop("gtf must contain transcript_id and gene_id")
  }
  
  gtf_df <- group_by(gtf_df, seqnames, strand, source, transcript_id, 
                     gene_id, gene_name)
  
  txs <- summarize(gtf_df, 
                   start = min(start),
                   end = max(end),
                   width = end - start + 1L,
                   type = "transcript",
                   score = NA, 
                   phase = NA,
                   exon_number = NA,
                   exon_id = NA)
  
  gtf_df <- group_by(gtf_df, seqnames, strand, source, gene_id, 
                     gene_name)
  
  genes <- summarize(gtf_df, 
                     start = min(start),
                     end = max(end),
                     width = end - start + 1L,
                     type = "gene",
                     score = NA, 
                     phase = NA,
                     transcript_id = NA,
                     exon_number = NA,
                     exon_id = NA)
  
  col_order <- quo(c("seqnames",
                     "start",
                     "end",
                     "width",
                     "strand",
                     "source",
                     "type",
                     "score",
                     "phase",
                     "gene_id",
                     "transcript_id",
                     "exon_number",
                     "exon_id",
                     "gene_name"))
  
  txs <- dplyr::select(txs, !!col_order)
  genes <- dplyr::select(genes, !!col_order)
  
  res <- suppressWarnings(bind_rows(list(gtf_df,
                                         txs,
                                         genes)))
  
  res <- arrange(res, seqnames, start)
  
  res <- ungroup(res)
  res
}


#' convert matrix to tibble with rownames
#' 
#' 
tidy_matrix <- function(df, row_ids = NULL){
  as_tibble(df, rownames = row_ids)
}


#' Exrtact utr's from a gtf file
#' 
#' 
#' 
get_3utrs <- function(gtf_file,
                      gene_ids,
                      gene_col = "gene_id",
                      transcript_col = "transcript_id") {
  
  gtf <- tidy_gtf(gtf_file)
  
  gtf <- filter(gtf, gene_id %in% gene_ids)
  
  gtf <- dplyr::filter(gtf,
                       type %in% c("exon", "CDS"))
  
  cds_gtf <- filter(gtf, type == "CDS") %>%
    group_by(transcript_id) %>%
    summarize(cds_start = unique(ifelse(strand == "+",
                                        min(start),
                                        max(end))),
              cds_end = unique(ifelse(strand == "+",
                                      max(end),
                                      min(start)))) %>%
    ungroup() %>% 
    left_join(gtf, ., by = "transcript_id")
  
  utr3_gtf <- cds_gtf %>%
    filter(type == "exon") %>%
    group_by(transcript_id) %>%
    mutate(utr3_exon = ifelse(
      strand == "+",
      ifelse(
        start >= cds_end,
        "yes",
        ifelse(cds_end >= start &
                 cds_end <= end,
               "overlapping",
               "no")
      ),
      ifelse(
        end <= cds_end,
        "yes",
        ifelse(cds_end >= start &
                 cds_end <= end,
               "overlapping",
               "no")
      )
    )) %>%
    ungroup() %>%
    filter(utr3_exon != "no")
  
  utr3_res <- utr3_gtf %>%
    group_by(transcript_id) %>%
    mutate(
      start = ifelse(strand == "+" & utr3_exon == "overlapping",
                     cds_end,
                     start),
      end = ifelse(strand == "-" & utr3_exon == "overlapping",
                   cds_end,
                   end)
    )  %>%
    filter(start != end) %>%
    select(chrom:transcript_id)
  
  utr3_res <- ungroup(utr3_res)
  utr3_res <- mutate(utr3_res, start = start - 1)
  utr3_res <- mutate_if(utr3_res, is.factor, as.character)
  res <- group_by(utr3_res, gene_id, strand) %>%
    bed_merge() %>%
    ungroup()
  
  res
}


get_tss <- function(gtf_file,
                    gene_ids = NULL){
  
  gtf <- tidy_gtf(gtf_file, zero_based_coords = TRUE)
  
  if(!is.null(gene_ids)){
    gtf <- filter(gtf,  gene_id %in% gene_ids)
  }
  
  tss <- filter(gtf, type == "gene") %>%
    mutate(end = ifelse(strand == "+",
                        start + 1,
                        end),
           start = ifelse(strand == "-",
                          end - 1,
                          start)) %>% 
    select(chrom, start, end, gene_name, gene_id, strand)
  
  tss
  
}
#' write out bed file without mangling numeric characters.
#' 
write_bed <- function(df, path = ""){
  res <- mutate_if(df, is.numeric, as.character)
  write_tsv(res, path, col_names = F)
}


#' add pre-mrna designation
add_pre <- function(chr_vec, pattern = "pre_"){
  str_c(pattern, chr_vec)
}  

#' remove pre-mrna designation 
no_pre <- function(chr_vec, pattern = "^pre_"){
  str_remove(chr_vec, pattern)
}  

#' calc zscore
zscore <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)}

#' 
plot_expr <- function(genes, 
                      cts, 
                      lib_to_sample, 
                      group = NULL,
                      facet = NULL){
  
  long_tidy <- cts[genes, , drop = FALSE] %>%
    as.data.frame(.) %>%
    rownames_to_column("gene") %>%
    gather(sample, expr, -gene)
  
  join_col <- colnames(lib_to_sample)[1]
  plt_dat <- left_join(long_tidy,
                       lib_to_sample, 
                       by = c("sample" = join_col))
  
  add_facet <- !is.null(facet)
  add_group <- !is.null(group)
  
  if(add_group && !add_facet){
    grp <- sym(group)
    plt_dat <- group_by(plt_dat, !!grp) %>% 
      summarize(expr = mean(expr)) %>% 
      ungroup() %>% 
      dplyr::rename(sample = !!grp)
  } else if (!add_group && add_facet){
    grp <- sym(group)
    fac <- sym(facet)
    plt_dat <- group_by(plt_dat, !!grp, !!fac) %>% 
      summarize(expr = mean(expr)) %>% 
      ungroup() %>% 
      dplyr::rename(sample = !!grp)
  }
  
  p <- ggplot(plt_dat, aes(sample, expr)) +
    geom_point() +
    labs(x = "",
         y = "Expression")
  
  if (add_facet){
    p <- p + facet_wrap(as.formula(paste("~", facet)))
  }
  
  if(length(genes) == 1){
    p <- p + labs(title = genes)
  }
  
  p + 
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}


plot_exon_intron_expr <- function(gene, 
                                  count_matrices,
                                  line_col = "blue",
                                  line_width = 2){
  if(length(gene) > 1){
    stop("only 1 gene at a time")
  }
  
  if(is.null(names(count_matrices))){
    names(count_matrices) <- c("exon", "intron")
  }
  long_tidy <- map_dfr(count_matrices, 
                       ~.x[gene, , drop = FALSE] %>%
                         as.data.frame(.) %>%
                         rownames_to_column("gene") %>%
                         gather(sample, expr, -gene),
                       .id = "count_type")
  
  p <- ggplot(long_tidy, aes(sample, expr, group = 1)) +
    geom_line(color = line_col, size = line_width) + 
    geom_point(color = "black") +
    ylim(c(0, NA)) + 
    facet_wrap(~count_type, scales = "free_y", ncol = 1) +
    labs(x = "",
         y = "TPM")
  
  
  p <- p + labs(title = gene)
  
  p + 
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}



#' emulate python zip
zip <- function(...) { Map(c, ...) }

get_bw <- function(bw_fn,
                   chrom = NULL,
                   start = NULL,
                   end = NULL,
                   bigWigToBedGraph = "/usr/local/bin/bigWigToBedGraph") {
  
  ucsc_args <- vector("character")
  
  if (!is.null(chrom)) {
    ucsc_args <- c(ucsc_args, paste0("-chrom=", chrom))
  }
  
  if (!is.null(start)) {
    ucsc_args <- c(ucsc_args, paste0("-start=", start))
  }
  
  if (!is.null(end)) {
    ucsc_args <- c(ucsc_args, paste0("-end=", end))
  }
  
  ucsc_args <- c(ucsc_args,
                 bw_fn,
                 "stdout")
  dat <- system2(bigWigToBedGraph,
                 ucsc_args,
                 stdout = TRUE)
  
  if(length(dat) > 0)
    res <- read_tsv(
      I(dat),
      col_names = c("contig",
                    "start",
                    "end",
                    "coverage")
    ) else{
      res <- tibble(contig = character(), 
                    start = numeric(), 
                    end = numeric(),
                    coverage = numeric())
    }
  
  res
}

#' Plot tx coverage
#' 
#' @param bams
#' @param features
#' @param gtf_obj

plot_tx_coverage <- function(bw_fns, 
                             transcripts,
                             gtf_obj,
                             scale_y = FALSE,
                             gnome = "dm6",
                             annotation_label = NULL,
                             ...){
  library(GenomicFeatures)
  library(Gviz)
  library(GenomicAlignments)
  
  options(ucscChromosomeNames=FALSE)
  
  #subset gtf to get min max coords to plot
  gfeatures <- gtf_to_gviz(gtf_obj, transcripts)
  
  chrom <- as.vector(unique(gfeatures$chromosome))[1]
  start <- min(gfeatures$start)
  end <- max(gfeatures$end)
  strand <- as.vector(unique(gfeatures$strand))[1]
  
  # get coverage
  cov_list <- map(bw_fns, 
                  ~get_bw(.x,
                          chrom = chrom,
                          start = start,
                          end = end) %>% 
                    makeGRangesFromDataFrame(., 
                                             keep.extra.columns = TRUE,
                                             seqnames.field = 'contig',
                                             starts.in.df.are.0based = TRUE,
                                             ignore.strand = TRUE))
  
  if(scale_y){
    ymax <- map(cov_list, ~.x$coverage) %>% unlist() %>% max()
    y_limits <- c(0, ymax)
  } else {
    y_limits <- NULL
  }
  
  lib_names <- names(bw_fns)
  
  if(is.null(lib_names)){
    lib_names <- as.character(1:length(bw_fns))
  } 
  names(cov_list) <- lib_names 
  
  dtracks <- pmap(list(cov_list,
                       names(cov_list),
                       viridis(length(cov_list))),
                  function(x, y, z){
                    DataTrack(range = x, 
                              name =  as.character(y),
                              chromosome = chrom, 
                              genome = gnome,
                              col.line = z,
                              ylim = y_limits
                    )})
  
  grtrack <- GeneRegionTrack(gfeatures,
                             genome = gnome, 
                             chromosome=chrom, 
                             name="",
                             #   stacking = "dense",
                             col.title = "black")
  options(scipen=16)
  if(is.null(annotation_label)){
    annotation_label <- unique(gfeatures$symbol)[1]
  }
  plotTracks(c(list(grtrack), 
               dtracks),
             chromosome = chrom,
             from = start,
             to = end,
             type = "S",
             reverseStrand = strand == "-", 
             col.title = "black",
             col.axis = "black",
             background.title ="transparent",
             lwd = 1.0, 
             innerMargin = 10,
             margin = 50,
             main = annotation_label,
             ...)
  
}

gtf_to_gviz <- function(gtf_obj, features){
  tx_granges <- gtf_obj[gtf_obj$transcript_id %in% features] 
  tx_df <- tx_granges %>% 
    as.data.frame() %>% 
    filter(type == "exon") %>% 
    dplyr::select(chromosome = seqnames,
                  start:strand,
                  gene = gene_id,
                  transcript = transcript_id,
                  exon = exon_id,
                  symbol = gene_name)
  tx_df
  
}


#' Set defaults in case of namespace conflicts
#' 
select <- function(...) {
  dplyr::select(...)
}

slice <- function(...){
  dplyr::slice(...)
}


#' Find exact non-overlapping matches in DNA seq
#' 
#' @param seq character vector of sequences to query  
#' @param q_str string to find (case insensitive query)
#' @param rev_comp also seach reverse complement (Default = FALSE)
#' @examples 
#' seq <- "ATCGTAGCTGATGAATAAATGCTGATGCTGTAG"
#' find_seq_match(seq, "AATAAA")

find_seq_match <- function(seq, q_str, rev_comp = FALSE) {
  seq <- str_to_upper(seq)
  q_str <- str_to_upper(q_str)
  stringr::str_locate_all(seq, fixed(q_str))
}

count_seq_matches <- function(seq, q_str, rev_comp = FALSE){
  map_int(find_seq_match(seq, q_str), nrow)
}

ggplot2::theme_set(theme_cowplot())


# calculate effective lengths for TPM measurements using kallisto algorithm
# salmon uses the same approach 
# rewritten in R based on kallisto c++ code
k_trunc_guassian_fld <- function(start, stop, mean_val, sd_val) {
  
  n = stop - start;
  mean_fl = vector("double", n);
  
  total_mass = 0.0
  total_density = 0.0
  
  lens <- 1:stop
  
  # z-scores 
  x <- (lens - mean_val) / sd_val
  
  # density 
  cur_density = exp(- 0.5 * x * x) / sd_val
  
  cumulative_den = cumsum(cur_density)
  cumultative_mass = cumsum(cur_density * lens)
  
  is_valid <- cumultative_mass > 0
  mean_fl[is_valid] <- cumultative_mass[is_valid] / cumulative_den[is_valid]
  
  names(mean_fl) <- lens
  mean_fl
}

k_get_frag_len_means <- function(tx_lens, frag_dis, max_frag_len = 800) {
  # use mean provided if larger than max frag length
  marginal_mean <- frag_dis[length(frag_dis)]
  tx_eff_len_means <- ifelse(tx_lens >= max_frag_len,
                             marginal_mean,
                             frag_dis[tx_lens])
  tx_eff_len_means
}


get_effective_lengths <- function(tx_lens, 
                                  mean_frag_len,
                                  sd_frag_len, 
                                  max_fragment_len = 800) {
  
  fld <- k_trunc_guassian_fld(0, max_fragment_len, mean_frag_len, sd_frag_len)
  tx_eff_lens_mean <- k_get_frag_len_means(tx_lens, fld )
  
  tx_eff_lens <- tx_lens - tx_eff_lens_mean + 1
  too_short <- tx_eff_lens < 1
  tx_eff_lens[too_short] <- tx_lens[too_short]
  
  tx_eff_lens
}

# example usage
# tx_lens <- sample(1:10000, 50000, replace = TRUE)
# 
# fld_mean <- 200
# fld_sd <- 20
# 
# get_effective_lengths(tx_lens, 
#                       fld_mean,
#                       fld_sd) %>% 
#   data.frame(tx_lens, .) %>% 
#   arrange(tx_lens) 


#'@export
calc_features <- function(gtf, features, gtf_field = "transcript_id"){
  
  if (typeof(gtf) == "character"){
    if (file.exists(gtf)){
      message("reading gtf file from ", gtf)
      gtf <- tidy_gtf(gtf)
    } else {
      stop("gtf file not found")
    }
  }
  
  if(!is_tidy_gtf(gtf)){
    stop("gtf input is not a valid tidy_gtf")
  }
  
  gtf_tbl <- filter_gtf(gtf, features, col = gtf_field)
  
  three_prime_utrs <- get_3utrs(gtf_tbl)
  five_prime_utrs <- get_5utrs(gtf_tbl)
  cds <- get_cds(gtf_tbl)
  res <- dplyr::bind_rows(three_prime_utrs, five_prime_utrs, cds)
  
  res <- gtf_to_bed6(res, score = "transcript_id")
  
  res
}

num_exons <- function(gtf){
  dplyr::group_by(gtf, transcript_id) %>%
    dplyr::summarize(n_exons = dplyr::n())
}

num_introns <- function(gtf){
  res <- num_exons(gtf)
  res <- dplyr::mutate(res, n_introns = n_exons - 1)
  dplyr::select(res, -n_exons)
}


tx_len <- function(gtf) {
  dplyr::mutate(gtf,
                flen = end - start) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarize(flen = sum(flen))
}

#'@export
calc_other_features <- function(gtf, features, gtf_field = "transcript_id") {
  if (typeof(gtf) == "character"){
    if (file.exists(gtf)){
      message("reading gtf file from ", gtf)
      gtf <- tidy_gtf(gtf)
    } else {
      stop("gtf file not found")
    }
  }
  
  if(!is_tidy_gtf(gtf)){
    stop("gtf input is not a valid tidy_gtf")
  }
  
  gtf_tbl <- filter_gtf(gtf, features, col = gtf_field)
  exon_tbl <- dplyr::filter(gtf_tbl, type == "exon")
  tx_tbl <- dplyr::filter(gtf_tbl, type == "transcript")
  tx_lengths <- tx_len(tx_tbl)
  n_exon <- num_exons(exon_tbl)
  n_intron <- num_introns(exon_tbl)
  res <- Reduce(function(x, y) dplyr::left_join(x,y, by = c("transcript_id")),
                list(tx_lengths, n_exon, n_intron))
  res
}

#' Calculate introns from a gtf file
#'
#'
#'
get_introns <- function(gtf) {
  
  gtf <- dplyr::filter(gtf,
                       type %in% c("exon"))
  res <- dplyr::group_by(gtf, transcript_id) %>%
    dplyr::arrange(start, .by_group = TRUE)
  
  res <- dplyr::mutate(res,
                       .start = end,
                       .end = lead(start),
                       score = ifelse(score == 1,
                                      1,
                                      score - 1),
                       start = .start,
                       end = .end) %>%
    dplyr::ungroup()
  
  res <- dplyr::select(res, -.start, -.end)
  res <- dplyr::filter(res,
                       !is.na(start),
                       !is.na(end),
                       start < end)
  
  res <- dplyr::arrange(res, chrom, start, end)
  
  res
}


get_tss <- function(gtf) {
  res <- dplyr::filter(gtf, type == "transcript")
  res <- mutate(res,
                start = ifelse(strand == "+",
                               start,
                               end - 1),
                end = ifelse(strand == "+",
                             start + 1,
                             end))
  # res <- select(res, chrom, start, end, strand)
  res
}


get_tts <- function(gtf, utr_id = "five_prime_utr") {
  res <- dplyr::filter(gtf, type == "transcript")
  res <- mutate(res,
                start = ifelse(strand == "+",
                               end - 1,
                               start),
                end = ifelse(strand == "+",
                             end,
                             start + 1))
  # res <- select(res, chrom, start, end, strand)
  res
}

#' Extract 3' utr's from a gtf file
#'
#'
#'
get_3utrs <- function(gtf, utr_id = "three_prime_utr") {
  
  if (any(utr_id %in% gtf$type)){
    res <- dplyr::filter(gtf, type == utr_id)
    return(res)
  }
  
  gtf <- dplyr::filter(gtf,
                       type %in% c("exon", "CDS"))
  
  cds_gtf <- filter(gtf, type == "CDS") %>%
    group_by(transcript_id) %>%
    summarize(cds_start = unique(ifelse(strand == "+",
                                        min(start),
                                        max(end))),
              cds_end = unique(ifelse(strand == "+",
                                      max(end),
                                      min(start)))) %>%
    ungroup() %>%
    left_join(gtf, ., by = "transcript_id")
  
  utr3_gtf <- cds_gtf %>%
    filter(type == "exon") %>%
    group_by(transcript_id) %>%
    mutate(utr3_exon = ifelse(
      strand == "+",
      ifelse(
        start >= cds_end,
        "yes",
        ifelse(cds_end >= start &
                 cds_end <= end,
               "overlapping",
               "no")
      ),
      ifelse(
        end <= cds_end,
        "yes",
        ifelse(cds_end >= start &
                 cds_end <= end,
               "overlapping",
               "no")
      )
    )) %>%
    ungroup() %>%
    filter(utr3_exon != "no")
  
  utr3_res <- utr3_gtf %>%
    group_by(transcript_id) %>%
    mutate(
      start = ifelse(strand == "+" & utr3_exon == "overlapping",
                     cds_end,
                     start),
      end = ifelse(strand == "-" & utr3_exon == "overlapping",
                   cds_end,
                   end)
    )  %>%
    filter(start != end) %>%
    select(chrom:transcript_id)
  
  utr3_res <- ungroup(utr3_res)
  utr3_res <- mutate(utr3_res, start = start - 1)
  utr3_res <- mutate_if(utr3_res, is.factor, as.character)
  res <- group_by(utr3_res, gene_id, strand) %>%
    bed_merge() %>%
    ungroup()
  
  res
}

#' Extract 3' utr's from a gtf file
#'
#'
#'
get_5utrs <- function(gtf, utr_id = "five_prime_utr") {
  
  if (any(utr_id %in% gtf$type)){
    res <- dplyr::filter(gtf, type == utr_id)
    return(res)
  }
  
}


#' Extract cds from a gtf file
#'
get_cds <- function(gtf, cds_id = "CDS") {
  
  if (any(cds_id %in% gtf$type)){
    res <- dplyr::filter(gtf, type == cds_id)
    return(res)
  }
  
}


#'@export
gtf_to_bed6 <- function(gtf,
                        name_col = "type",
                        score_col = "exon_number") {
  cols <- c(
    "chrom",
    "start",
    "end",
    name_col,
    score_col,
    "strand"
  )
  
  gtf[, cols]
}

filter_gtf <- function(gtf, features, col = "transcript_id"){
  if(!is_tidy_gtf(gtf)){
    stop("gtf input incorrect")
  }
  
  if (!all(features %in% gtf[[col]])) {
    stop("unable to find some features")
  } else {
    res <- dplyr::filter(gtf, !!rlang::sym(col) %in% features)
  }
  
  return(res)
}

#'@export
is_tidy_gtf <- function(x){
  
  cols <- c("chrom",
            "start",
            "end",
            "width",
            "strand",
            "source",
            "type",
            "score",
            "phase",
            "gene_id",
            "transcript_id",
            "exon_number",
            "exon_id")
  
  right_cols <- all(cols %in% colnames(x))
  
  right_cols
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


#'Extract sequences from fasta using gtf input
#'
#'@param gtf dataframe containing gtf style entries, either in as bed (chrom, start, end) or (seqnames)
#'@param fasta_path path to fasta file with .fai index
#'@return data_frame containing transcript id column and sequence column
#'@export
gtf_to_seq <- function(gtf, fasta_path, type_field = "exon"){
  
  df_cols <- c("chrom", "start", "end",
               "transcript_id", "exon_number", "strand")
  
  gtf_cols <- c("seqnames", "start", "end",
                "transcript_id", "exon_number", "strand")
  
  if (!all(df_cols %in% colnames(gtf))) {
    if (all(gtf_cols %in% colnames(gtf))){
      df <- dplyr::rename(df, chrom = seqnames)
      df <- dplyr::mutate(df, start = start - 1)
    } else {
      stop("unknown columns in supplied dataframe")
    }
  }
  
  df <- dplyr::filter(gtf, type == type_field)
  df <- dplyr::select(df,
                      chrom,
                      start,
                      end,
                      transcript_id,
                      exon_number,
                      strand)
  
  if(!"kentr" %in% installed.packages()){
    stop("please install the kentr package",
         "  devtools::install_github('kriemo/kentr')")
  }
  
  seq_df <- kentr::get_sequences(df, fasta_path, strand = FALSE)
  seq_df <- dplyr::group_by(seq_df, transcript_id)
  seq_df <- dplyr::arrange(seq_df, chrom, start, end, .by_group = TRUE)
  res <- dplyr::summarize(seq_df,
                          strand = unique(strand),
                          seq = stringr::str_c(seq, collapse = ""))
  res <- dplyr::mutate(res,
                       seq = ifelse(strand == "+",
                                    seq,
                                    kentr::revComp(seq)))
  res
}

#' From biostrings
#'@export
GENETIC_CODE <- c(
  TTT="F",
  TTC="F",
  TTA="L",
  TTG="L",
  
  TCT="S",
  TCC="S",
  TCA="S",
  TCG="S",
  
  TAT="Y",
  TAC="Y",
  TAA="*",
  TAG="*",
  
  TGT="C",
  TGC="C",
  TGA="*",
  TGG="W",
  
  CTT="L",
  CTC="L",
  CTA="L",
  CTG="L",
  
  CCT="P",
  CCC="P",
  CCA="P",
  CCG="P",
  
  CAT="H",
  CAC="H",
  CAA="Q",
  CAG="Q",
  
  CGT="R",
  CGC="R",
  CGA="R",
  CGG="R",
  
  ATT="I",
  ATC="I",
  ATA="I",
  ATG="M",
  
  ACT="T",
  ACC="T",
  ACA="T",
  ACG="T",
  
  AAT="N",
  AAC="N",
  AAA="K",
  AAG="K",
  
  AGT="S",
  AGC="S",
  AGA="R",
  AGG="R",
  
  GTT="V",
  GTC="V",
  GTA="V",
  GTG="V",
  
  GCT="A",
  GCC="A",
  GCA="A",
  GCG="A",
  
  GAT="D",
  GAC="D",
  GAA="E",
  GAG="E",
  
  GGT="G",
  GGC="G",
  GGA="G",
  GGG="G"
)

#'@export
translate_seq <- function(vec){
  stringr::str_match_all(vec, ".{3}") %>%
    purrr::map_chr(~str_c(GENETIC_CODE[.x], collapse = ""))
}


is_na_equal <- function(x, y, ignore_na = TRUE){
  x_na <- is.na(x)
  y_na <- is.na(y)
  
  if(!any(x_na, y_na)){
    return(x == y)
  } else if(x_na && y_na){
    return(NA)
  } else {
    return(FALSE)
  }
  
}

comp_feature_seq <- function(x, y, gtf_df, ftype, fasta_path) {
  ref <- tibble(
    transcript_id = as.vector(rbind(x, y)),
    tx_type = rep(c("query_tx", "ref_tx"),
                  length(x)),
    id = rep(1:length(x), each = 2)
  )
  
  feats <- filter(gtf_df,
                  transcript_id %in% ref$transcript_id,
                  type == ftype) %>%
   gtf_to_seq(.,
                          fasta_path = fasta_path,
                          type_field = ftype) %>%
    dplyr::select(-strand)
  feats <- left_join(ref, feats, by = c("transcript_id"))
  
  out_bool <- sym(paste0("same_", ftype))
  out_len <- sym(paste0(ftype, "_length"))
  
  feats <- group_by(feats, id) %>%
    mutate(!!out_bool := is_na_equal(dplyr::first(seq),dplyr::nth(seq, 2))) %>%
    ungroup() %>%
    mutate(!!out_len := nchar(seq)) %>%
    dplyr::select(-seq,-id)
  
  feats
}

#'@export
comp_feature_seqs <- function(x, y,
                              gtf_df,
                              ftypes = c("three_prime_utr",
                                         "CDS",
                                         "five_prime_utr"),
                              fasta_path = "~/Projects/shared_dbases/genomes/drosophila/Drosophila_melanogaster.BDGP6.dna.toplevel.fa"){
  
  if(length(x) != length(y)){
    stop("only equal length x and y transcript vectors can be compared")
  }
  
  tx_gtf <- dplyr::filter(gtf_df, transcript_id %in% c(x,y))
  
  names(ftypes) <- ftypes
  res <- map(ftypes, ~comp_feature_seq(x, y, tx_gtf, .x, fasta_path))
  res <- Reduce(function(x, y) dplyr::inner_join(x, y,
                                                 by = c("transcript_id",
                                                        "tx_type")),
                res)
  res <- dplyr::select(res,
                       transcript_id,
                       tx_type,
                       starts_with("same_"),
                       ends_with("_length"))
  res
}



comp_feature_pos <- function(x, y,
                             gtf_df,
                             ftype){
  
  ref <- tibble(
    transcript_id = as.vector(rbind(x, y)),
    tx_type = rep(c("query_tx", "ref_tx"),
                  length(x)),
    id = rep(1:length(x), each = 2)
  )
  
  feats <- filter(gtf_df,
                  transcript_id %in% ref$transcript_id)
  
  get_feat <- switch(ftype,
                     "introns" = get_introns,
                     "tss" = get_tss,
                     "tts" = get_tts)
  
  feats <- get_feat(feats)
  
  feats <- left_join(ref, feats, by = c("transcript_id"))
  
  out_bool <- sym(paste0("same_", ftype))
  
  feats_same <- dplyr::select(feats, tx_type, id, start, end) %>%
    split(., .$id)  %>%
    purrr::map_df(function(x) {
      split_vec <- split(x, x$tx_type)
      same_start <- identical(sort(split_vec[[1]]$start), sort(split_vec[[2]]$start))
      same_end <- identical(sort(split_vec[[1]]$end), sort(split_vec[[2]]$end))
      tibble(!!out_bool := all(c(same_start, same_end)))},
      .id = "id") %>%
    mutate(id = as.integer(id))
  
  feats <- left_join(feats, feats_same, by = "id")
  feats
}




#'@export
comp_feature_positions <- function(x, y,
                                   gtf_df,
                                   ftypes = c("tss",
                                              "introns",
                                              "tts")){
  
  if(length(x) != length(y)){
    stop("only equal length x and y transcript vectors can be compared")
  }
  
  tx_gtf <- dplyr::filter(gtf_df, transcript_id %in% c(x,y))
  
  names(ftypes) <- ftypes
  res <- map(ftypes, ~comp_feature_pos(x, y, tx_gtf, .x))
  res <- Reduce(function(x, y) dplyr::inner_join(x, y,
                                                 by = c("transcript_id",
                                                        "tx_type")),
                res)
  res <- dplyr::select(res,
                       transcript_id,
                       tx_type,
                       starts_with("same_"),
                       ends_with("_length")) %>%
    unique()
  res
}





#### Annotations ####

annots <- list(
  drosophila = list(
    gtf  = file.path(db_dir, "drosophila", "ext", "Drosophila_melanogaster.BDGP6.84.gtf"),
    genome = file.path(db_dir, "drosophila", "ext", "Drosophila_melanogaster.BDGP6.dna.toplevel.fa"),
    chroms = file.path(db_dir, "drosophila", "ext", "Drosophila_melanogaster.BDGP6.dna.toplevel.fa.fai")
  ),
  coral = list(
    gtf  = file.path(db_dir, "coral", "ext", "Montipora_capitata_HIv2.genes.gtf"),
    genome = file.path(db_dir, "coral", "ext", "Montipora_capitata_HIv2.assembly.fasta"),
    chroms = file.path(db_dir, "coral", "ext", "Montipora_capitata_HIv2.assembly.fasta.fai")
  ),
  xenopus_t = list(
    gtf  = file.path(db_dir, "xenopus_t", "ext", "xenTro9.gtf"),
    genome = file.path(db_dir, "xenopus_t", "ext", "xenTro9.fa"),
    chroms = file.path(db_dir, "xenopus_t", "ext", "xenTro9.fa.fai")
  ),
  zebrafish = list(
    gtf  = file.path(db_dir, "zebrafish", "ext", "danRer10.gtf"),
    genome = file.path(db_dir, "zebrafish", "ext", "danRer10.fa"),
    chroms = file.path(db_dir, "zebrafish", "ext", "danRer10.fa.fai")
  )
)












