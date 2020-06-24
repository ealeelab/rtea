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
  make_option(c("-b", "--bamfile"), help = "bam file path"),
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

library(Rsamtools)
chrnames <- names(scanBamHeader(bamfile)[[1]]$targets)
rtea <- fread(rteafile)
if(startsWith(chrnames[1], "chr")) {
  rtea[, chr := pastechr(chr)]
} else {
  rtea[, chr := sub("^chr", "", chr)]
}
cnt <- countClippedReads.ctea(rtea[, chr:TEbreak], bamfile, 
                              strandedness = "non-stranded",
                              fusiontype_cutoff = Inf)
cnt %<>% data.table(rtea[, !names(.), with = F])
fwrite(cnt, outfile, sep = "\t", na = "NA", quote = F)
