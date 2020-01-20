library(data.table)
library(magrittr)
require(Rcpp)

thisdir <- getwd()
if(!file.exists(file.path(thisdir, "ctea_functions.R"))) {
  warning("It is recommended to parse this file using sys.source with chdir = T")
}

bamtoolsinclude <- "/home/boramlee/tools/bamtools-2.5.1/include/bamtools"
bamtoolslib <- "/data2/boramlee/tools/bamtools-2.5.1/lib64"
clipcpp <- file.path(thisdir, "clippedSequence.cpp")
combinecpp <- file.path(thisdir, "combineNeighbor.cpp")
awkfile <- file.path(thisdir, "teaify-0.7-no-combine.awk")
if(!exists("refFa"))
  refFa <- file.path(thisdir, "../ref/ctea/repeat_LINE1_ALU_SVA_HERV_human_youngTE.fa")
if(!grepl("bamtools", Sys.getenv("PKG_CXXFLAGS"))) 
  Sys.setenv(PKG_CXXFLAGS = paste0("-I", bamtoolsinclude))
if(!grepl("bamtools", Sys.getenv("PKG_LIBS"))) 
  Sys.setenv(PKG_LIBS = paste0("-L", bamtoolslib, " -lbamtools"))
dyn.load(file.path(bamtoolslib, "libbamtools.so"))
sourceCpp(clipcpp)
sourceCpp(combinecpp)


runctea <- function(bamfile, outprefix, refFa = refFa, threads = 2, removetmp = T) {
  setDTthreads(threads)
  tmpdir <- paste0(outprefix, ".tmp")
  dir.create(tmpdir, recursive = T)
  clipprefix <- file.path(tmpdir, "01-clipped")
  clipped_f <- paste0(clipprefix, "-f-sorted")
  clipped_r <- paste0(clipprefix, "-r-sorted")
  contigs_f <- file.path(tmpdir, "02-clipped-f.fa")
  contigs_r <- file.path(tmpdir, "02-clipped-r.fa")
  sam_f <- file.path(tmpdir, "03-clipped-f.sam")
  sam_r <- file.path(tmpdir, "03-clipped-r.sam")
  filtered_f <- file.path(tmpdir, "04-filtered-f")
  filtered_r <- file.path(tmpdir, "04-filtered-r")
  sorted_f <- file.path(tmpdir, "04-filtered-f-sorted")
  sorted_r <- file.path(tmpdir, "04-filtered-r-sorted")
  combined_f <- file.path(tmpdir, "05-combined-f")
  combined_r <- file.path(tmpdir, "05-combined-r")
  cteafile <- paste0(outprefix, ".ctea")
  
  clip(bamfile, clipprefix)
  contigen(clipped_f, contigs_f)
  contigen(clipped_r, contigs_r)
  bwa_aln(refFa, contigs_f, sam_f, threads = threads)
  bwa_aln(refFa, contigs_r, sam_r, threads = threads)
  filter_family(sam_f, filtered_f)
  filter_family(sam_r, filtered_r)
  combine(sorted_f, combined_f, "f")
  combine(sorted_r, combined_r, "r")
  give_refName(combined_f, combined_r, cteafile, bamfile)
  
  if(removetmp) {
    unlink(tmpdir, recursive = T)
  }
}

contigen <- function(infile, outfile) {
  clipseq <- fread(infile)
  setnames(clipseq, c("refID", "pos", "seq"))
  writeLines(clipseq[, paste0(">", refID, ";", pos, ";", seq, "\n", seq)], outfile)
}

bwa_aln <- function(refFa, infile, outfile, threads = 2) {
  saifile <- sub("[.]sam$", ".sai", outfile)
  cmd1 <- paste("bwa aln",
               "-t", threads,
               refFa,
               infile, ">",
               saifile)
  cmd2 <- paste("bwa samse -n 100",
                refFa,
                saifile,
                infile, ">",
                outfile)
  writeLines(cmd1)
  system(cmd1)
  writeLines(cmd2)
  system(cmd2)
}

filter_family <- function(infile, outprefix) {
  cmd1 <- paste("awk -f", awkfile,
                infile, ">",
                outprefix)
  cmd2 <- paste("sort -g -k1 -k2",
                outprefix, ">",
                paste0(outprefix, "-sorted"))
  writeLines("Filtering reads without family...")
  system(cmd1)
  writeLines(paste("Sorting", outprefix))
  system(cmd2)
}

give_refName <- function(ffile, rfile, outfile, bamfile) {
  refName <- get_refName(bamfile)
  ctea <- rbind(fread(ffile), fread(rfile))
  ctea[, refID := factor(refID, levels = refName, labels = names(refName))]
  setnames(ctea, "refID", "chr")
  fwrite(ctea, outfile, sep = "\t", quote = F)
}
