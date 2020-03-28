#!/bin/env Rscript

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
  make_option(c("-r", "--rteafile"), help = "rtea output file path"),
  make_option(c("-b", "--bamfile"), help = "wgs bam file path"),
  make_option(c("-o", "--outfile"), help = "output file path"),
  make_option(c("-t", "--threads"), type = "integer", default = 2, help = "number of threads"),
  make_option(c("--rtea_script"), default = file.path(thisdir, "../rtea_functions.R"), help = "rtea Rscript file")
)
opt <- parse_args(OptionParser(option_list = option_list))
attach(opt)

options(error = function(){
  sink(stderr())
  traceback(3)
  sink()
  save.image(paste0(outfile, ".err.RData"))
  q("no", status = 1, runLast = F)
})

options(mc.cores = threads)
setDTthreads(threads)

RTEA <- new.env()
sys.source(rtea_script, RTEA, chdir = T)
attach(RTEA)

rtea <- fread(rteafile)
rtea[, chr := paste0("chr", chr)]
wgscnt <- countClippedReads.ctea(rtea[, chr:TEbreak], bamfile)
wgscnt %<>% data.table(rtea[, gene_id:polyTE])
fwrite(wgscnt, outfile, sep = "\t", na = "NA", quote = F)
