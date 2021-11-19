#!/usr/bin/env Rscript
# Taken from https://gist.github.com/LTLA/6fde8edf558a0d32314043ab59e253cf 
# good example of how to interate through a file in R.
# reads assigned to both introns and exons will be ignored for our purposes

library(Rsubread)

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 5){
  stop("missing args: bam intron_saf exon_saf strand outfilename threads")
}

bam <- args[1]
isaf_fname <- args[2]
esaf_fname <- args[3]
strand <- args[4]
out_fname <- args[5]
threads <- args[6]

outdir <- dirname(out_fname)
sample_id <- gsub("_sorted.bam", "", basename(bam))

IX <- read.table(isaf_fname, header = TRUE, stringsAsFactors = FALSE)
EX <- read.table(esaf_fname, header = TRUE, stringsAsFactors = FALSE)

intron.counts <- featureCounts(bam, 
                               annot.ext=IX, 
                               reportReads="CORE", 
                               minOverlap=10,
                               reportReadsPath = outdir,
                               nthreads = threads,
                               strandSpecific = strand)

tmp_out <- file.path(outdir, paste0(basename(bam), ".featureCounts"))
intron_out <- file.path(outdir, paste0(sample_id, ".intron"))
file.rename(tmp_out, intron_out)

exon.counts <- featureCounts(bam, 
                             annot.ext=EX, 
                             reportReads="CORE", 
                             fracOverlap=1.0,
                             reportReadsPath = outdir,
                             nthreads = threads,
                             strandSpecific = strand)
exon_out <- file.path(outdir, paste0(sample_id, ".exon"))
file.rename(tmp_out, exon_out)

ugenes <- unique(c(unique(EX$GeneId), unique(IX$GeneID)))
out.exon <- out.intron <- out.both <- setNames(integer(length(ugenes)), ugenes)

fin <- file(intron_out, open='r')
fex <- file(exon_out, open='r')
n <- 0
n_lines <- 1e6
repeat {
    n <- n + 1e6
    curin <- read.table(fin, nrow=n_lines, colClasses=c("character", "character", "integer", "character"), stringsAsFactor=FALSE)
    curex <- read.table(fex, nrow=n_lines, colClasses=c("character", "character", "integer", "character"), stringsAsFactor=FALSE)

    # Figuring out what to do with reads based on their mutual locations.
    intron <- !is.na(curin[,4]) 
    exon <- !is.na(curex[,4])

    intron.only <- intron & !exon
    exon.only <- exon & !intron
    both <- intron & exon & curin[,4]==curex[,4]

    curin.counts <- table(curin[intron.only,4])
    curex.counts <- table(curex[exon.only,4])
    both.counts <- table(curin[both,4])

    out.intron[names(curin.counts)] <- out.intron[names(curin.counts)] + curin.counts
    out.exon[names(curex.counts)] <- out.exon[names(curex.counts)] + curex.counts
    out.both[names(both.counts)] <- out.both[names(both.counts)] + both.counts

    if (nrow(curin)<n_lines) {
        break
    }
}

close(fex)
close(fin)

res <- data.frame(gene_id = names(out.intron),
                  intron_counts = out.intron,
                  exon_counts = out.exon[names(out.intron)], 
                  row.names = NULL)

write.table(res, out_fname, sep = "\t", row.names = FALSE)
write.table(intron.counts$annotation, 
            file.path(outdir, 
                      paste0(sample_id, "_intron_annot_out.tsv")), 
            sep = "\t",
            row.names = FALSE)
write.table(exon.counts$annotation, 
            file.path(outdir, 
                      paste0(sample_id, "_exon_annot_out.tsv")), 
            sep = "\t",
            row.names = FALSE)
file.remove(exon_out)
file.remove(intron_out)   
