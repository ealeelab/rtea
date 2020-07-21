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
  make_option(c("-c", "--cteafile"), help = "ctea output file path"),
  make_option(c("-o", "--outfile"), help = "output directory"),
  make_option(c("-t", "--threads"), type = "integer", default = 2, help = "number of threads"),
  make_option(c("--genome_build"), default = "hg38", help = "reference genome build"),
  make_option(c("--rtea_script"), default = file.path(thisdir, "rtea_functions.R"), help = "rtea Rscript file"),
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

ctea <- readctea(cteafile)
ctea <- filterUnlocalized.ctea(ctea)
ctea <- filterSimpleRepeat.ctea(ctea)
ctea <- repeatPositon.ctea(ctea)
ctea <- filterSimpleSite(ctea)
ctea <- filterNoClip.ctea(ctea, threads = threads)
ctea <- TEcoordinate(ctea, threads = threads)
ctea[isPolyA == T, class := "PolyA"]
ctea <- ctea[isPolyA | TEscore > 0, ]
fwrite(ctea, outfile, sep = "\t", quote = F, na = "NA")
