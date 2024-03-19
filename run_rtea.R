#!/usr/bin/env Rscript

options(error = function(){
  sink(stderr())
  traceback(3)
  sink()
  save.image("cteaFilterErr.RData")
  q("no", status = 1, runLast = F)
})

library(magrittr)
library(data.table)
library(stringr)
library(optparse)

this <- commandArgs() %>%
  str_extract("(?<=file=).*") %>%
  na.omit
thisdir <- if(length(this) == 1) {
  dirname(this)
} else {
  getwd()
}

option_list <- list(
  make_option(c("-b", "--bamfile"), help = "bam file path"),
  make_option(c("-c", "--cteafile"), help = "ctea output file path"),
  make_option(c("-o", "--outfile"), help = "output file"),
  make_option(c("-t", "--threads"), type = "integer", default = 2, help = "number of threads"),
  make_option(c("--genome_build"), default = "hg38", help = "reference genome build"),
  make_option(c("--rtea_script"), default = file.path(thisdir, "rtea_functions.R"), help = "rtea Rscript file"),
  make_option(c("--filter_script"), default = file.path(thisdir, "ctea_filter.R"), help = "ctea_filter Rscript file"),
  make_option(c("--refdir"), help = "directory containing reference files")
)
opt <- parse_args(OptionParser(option_list = option_list))
attach(opt)
if(is.null(refdir)) refdir <- file.path(thisdir, "ref", genome_build)

options(
  mc.cores = threads,
  genome_build = genome_build,
  refdir = refdir
)
setDTthreads(threads)
RTEA <- new.env()
sys.source(rtea_script, RTEA, chdir = T)
attach(RTEA)

precntfile <- paste0(dirname(outfile), "/precount_", basename(outfile))
system2("Rscript",
        c(filter_script,
          "-c", cteafile,
          "-o", precntfile,
          "-t", threads,
          "--genome_build", genome_build,
          "--rtea_script", rtea_script,
          "--refdir", refdir
        )
)
ctea <- fread(precntfile)

rtea <- countClippedReads.ctea(ctea, bamfile, threads = threads)
rtea <- rtea[trueCnt >= 3 | (isPolyA & polyAcnt >=3), ]
rtea <- calculateDepth.rtea(rtea, bamfile)
rtea[, VAF := matchCnt / depth]
rtea <- annotate.ctea(rtea)
rtea <- polyTElocation.ctea(rtea)
rtea <- localHardClip(rtea, threads = threads)
rtea <- fusiontype(rtea)
rtea <- cntFilter.ctea(rtea)

colSelect <- c("chr", "pos", "ori", "class", "seq",
               "isPolyA", "posRepFamily", "posRep", 
               "TEfamily", "TEscore", "TEside", "TEbreak",
               "depth", "matchCnt", "VAF",
               "polyAcnt", "baseQual", "lowMapQual", "mateDist",
               "overhang", "gap", 
               "secondary", "nonspecificTE", "r1pstrand",
               "strand", "pos_type", "polyTE",
               "hardstart", "hardend", "hardTE", "hardDist",
               "fusion_tx_id", "tx_support_exon", "tx_support_intron",
               "fusion_type", "fusion_tx_biotype", "fusion_gene_id", "fusion_gene_name",
               "Filter")
rtea <- rtea[, colSelect, with = F]

writeLines(paste("Writing result to", outfile))
print(head(rtea))
fwrite(rtea, outfile, sep="\t", na="NA", quote=F)
