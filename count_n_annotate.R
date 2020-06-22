#!/bin/env Rscript

options(error = function(){
  sink(stderr())
  traceback(3)
  sink()
  save.image("rteaErr.RData")
  q("no", status = 1, runLast = F)
})

library(magrittr)
library(data.table)
library(parallel)
library(stringr)

load(commandArgs(T)[1])
if(!(exists("RTEA") & exists("opt") & exists("ctea"))) {
  stop("Input data is not from ctea_filter.R")
}
attach(RTEA)
attach(opt)

options(
  mc.cores = threads,
  genome_build = genome_build,
  refdir = refdir
)

rtea <- countClippedReads.ctea(ctea, bamfile, threads = threads)
rtea <- rtea[trueCnt >= 3 | (isPolyA & polyAcnt >=3), ]
rtea <- annotate.ctea(rtea)
rtea <- polyTElocation.ctea(rtea)
rtea <- localHardClip(rtea, threads = threads)
rtea <- fusiontype(rtea)
rtea <- cntFilter.ctea(rtea)

outfile <- paste0(outdir, "/", sample_name, ".rtea.txt")
writeLines(paste("Writing result to", outfile))
print(head(rtea))
fwrite(rtea, outfile, sep="\t", na="NA", quote=F)

