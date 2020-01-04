#!/usr/bin/Rscript
#SBATCH -p short
#SBATCH -t 0-01:30
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH -o rtea_pipeline.sh.%j.out
#SBATCH -e rtea_pipeline.sh.%j.err
#SBATCH --mail-type=FAIL

# 190829 filter if clipped sequence does not map to consensus sequence

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

option_list <- list( 
  make_option(c("-b", "--bam"), help = "bam file path"),
  make_option(c("-c", "--ctea"), help = "ctea output file path"),
  make_option(c("-s", "--scallop"), help = "scallop output file path"),
  make_option(c("-o", "--out"), help = "output file path"),
  make_option(c("-t", "--threads"), type = "integer", default = 2, help = "number of threads"),
  make_option(c("--build"), default = "hg19", help = "reference genome build"),
  make_option(c("--rtea"), default = file.path(thisdir, "rtea_functions.R"), help = "rtea Rscript file"),
  make_option(c("--refdir"), help = "directory containing reference files")
)
opt <- parse_args(OptionParser(option_list = option_list))
cteafile <- opt$ctea
bamfile <- opt$bam
scallopfile <- opt$scallop
outfile <- opt$out
threads <- opt$threads
genome_build <- opt$build
TEI_expression_Rscript <- opt$rtea
refdir <- opt$refdir
if(is.null(refdir)) refdir <- file.path(thisdir, "ref", genome_build)


options(
   threads = threads,
   genome_build = genome_build,
   refdir = refdir
)
sys.source(TEI_expression_Rscript, .GlobalEnv, chdir = T)

library(GenomicAlignments)
setDTthreads(threads)

ctea <- readctea(cteafile, threads = threads) %>%
        filterUnlocalized.ctea %>%
        filterSimpleRepeat.ctea %>%
        filterTEsite.ctea %>%
        filterNoClip.ctea(threads = threads) %>%
        TEcoordinate
ctea %<>% .[isPolyA | TEscore > 0, ]
ctea %<>% countClippedReads.ctea(bamfile, threads = threads)
rtea <- ctea[trueCnt >= 3 | (isPolyA & polyAcnt >=3), ]
rtea %<>% cntFilter.ctea %>%
          annotate.ctea %>%
          polyTElocation.ctea %>%
          localHardClip
if(!is.null(opt$scallop)) {
  rtea %<>% matchScallop.ctea(scallopfile)
}

writeLines(paste("Writing to", outfile))
saveRDS(rtea, "rtea.rds")
print(head(rtea))
fwrite(rtea, outfile, sep="\t", na="NA", quote=F)

