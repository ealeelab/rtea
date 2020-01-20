#!env Rscript

options(error = function(){
  sink(stderr())
  traceback(3)
  sink()
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

ctea_functions_Rscript <- file.path(thisdir, "ctea/ctea_functions.R")

option_list <- list( 
  make_option(c("-b", "--bam"), help = "bam file path"),
  make_option(c("-o", "--out"), help = "outdir/prefix"),
  make_option(c("-r", "--ref"), 
              default = file.path(thisdir, "ref/ctea/repeat_LINE1_ALU_SVA_HERV_human_youngTE.fa"),
              help = "TE reference fasta file"),
  make_option(c("-t", "--threads"), type = "integer", default = 2, help = "number of threads"),
  make_option(c("--keeptmp"), action = "store_true", default = FALSE, help = "do not remove tmp directory"),
  make_option(c("--include"), help = "bamtools include path"),
  make_option(c("--library"), help = "bamtools library path")
)
opt <- parse_args(OptionParser(option_list = option_list))
bamfile <- opt$bam
outprefix <- opt$out
refFa <- opt$ref
threads <- opt$threads
removetmp <- !opt$keeptmp
bamtoolsinclude <- opt$include
bamtoolslib <- opt$library

if(!is.null(bamtoolsinclude)) 
  Sys.setenv(PKG_CXXFLAGS = paste0("-I", bamtoolsinclude))
if(!is.null(bamtoolslib)) {
  Sys.setenv(PKG_LIBS = paste0("-L", bamtoolslib, " -lbamtools"))
  dyn.load(file.path(bamtoolslib, "libbamtools.so"))
}

CTEA <- new.env()
sys.source(ctea_functions_Rscript, CTEA, chdir = T)
attach(CTEA)
runctea(bamfile = bamfile, 
        outprefix = outprefix, 
        refFa = refFa, 
        threads = threads, 
        removetmp = removetmp)

