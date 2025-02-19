require(magrittr)
require(data.table)
library(parallel)
library(stringr)

thisdir <- getwd()
if(!file.exists(file.path(thisdir, "rtea_functions.R"))) {
  warning("It is recommended to parse this file using sys.source with chdir = T")
}

##### reference files #####
build <- getOption("genome_build", "hg38")
refdir <- getOption("refdir", file.path(thisdir, "ref", build))
polyTE_file <- file.path(refdir, "GRanges_union_non_reference_TEI.rds")
bsgenomePKG <- paste0("BSgenome.Hsapiens.UCSC.", build)
edbPKG <- switch(build,
                 hg38 = "EnsDb.Hsapiens.v86",
                 hg19 = "EnsDb.Hsapiens.v75"
)
gene_file <- file.path(refdir, "annotation_data.rds")
rmsk_file <- file.path(refdir, "rmsk.rds")
refTEfa <- getOption("refTEfa", file.path(refdir, "consensus_L1_ALU_SVA_HERV.fa"))

##### functions #####
msg <- function(...) {
  message(format(Sys.time(), "\n[%Y-%m-%d %X %Z] "), ...)
}

reattach <- function() {
  if(!exists("TEI")) TEI <- new.env()
  if("TEI" %in% search()) detach(TEI)
  sys.source(file.path(thisdir, "rtea_functions.R"), TEI, chdir = T)
  attach(TEI)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
  modev <- uniqv[which.max(tabulate(match(v, uniqv)))]
  attributes(modev) <- attributes(v)
  modev
}

pastechr <- function(x) {
  sub("^(chr)?", "chr", sub("MT", "M", x))
}

readctea <- function(cteafile, count_cutoff = 3, confidence_cutoff = 2) {
  # require(Biostrings)
  # require(parallel)
  msg(paste("Reading", cteafile))
  awkcmd <- paste("awk -F '\\t' '$4 >=", count_cutoff,
                  "&& $6 >=", confidence_cutoff,
                  "'",
                  cteafile
  )
  ctea <- fread(cmd = awkcmd, sep = "\t")
  setnames(ctea,
           c("chr", "pos", "ori", "cnt", "family", "confidence",
             "moreFamily", "seq", "morePos", "moreSeq")
           )
  # revseq <- ctea[["seq"]] %>%
  #   DNAStringSet %>%
  #   reverseComplement %>%
  #   as("character")
  # moreseq <- strsplit(ctea[["moreSeq"]], ",")
  # morepos <- strsplit(ctea[["morePos"]], ",")
  # seqmatch <- mcmapply(match, ctea[["seq"]], moreseq, mc.cores = threads)
  # revseqmatch <- mcmapply(match, revseq, moreseq, mc.cores = threads)
  # seqmatch[is.na(seqmatch)] <-revseqmatch[is.na(seqmatch)]
  # ctea[, pos := mcmapply(`[`, morepos, seqmatch, mc.cores = threads) %>% as.integer]
  # ctea[, seq := mcmapply(`[`, moreseq, seqmatch, mc.cores = threads)]
  ctea[, class := NA_character_]
  ctea[family %in% c("L1HS", "LINE1"), class := "L1"]
  ctea[toupper(substr(family, 1, 3)) == "ALU", class := "Alu"]
  ctea[substr(family, 1, 3) == "SVA", class := "SVA"]
  ctea[substr(family, 1, 3) %in% c("HER", "MER", "LTR", "MLT"), class := "HERV"]

  ctea[, .(chr, pos, ori, cnt, class, family, moreFamily, confidence, seq)]
}

readxtea <- function(xteafile, total_count_cutoff = 5) {
  msg(paste("Reading", xteafile))
  xtea <- fread(xteafile)
  colnm <- switch(as.character(ncol(xtea)),
                  "6" = c("chr", "pos", "lt", "rt", "mate", "class"),
                  "8" = c("chr", "pos", "lt", "rt", "mate", "dl", "dr", "class")
  )
  setnames(xtea, colnm)

  xtea[, total := pmax(lt + mate, rt + mate)]
  xtea <- xtea[total >= total_count_cutoff, ]
  xtea[, ori := "*"]
  xtea[lt == 0, ori := "r"]
  xtea[rt == 0, ori := "f"]
  data.table(xtea[, .(chr, pos, ori, total)], xtea[, !c("chr", "pos", "ori", "total")])
}

filterUnlocalized.ctea <- function(ctea) {
  ctea[sub("^chr", "", chr) %in% c(1:22, c("X", "Y")), ]
}

filterSimpleRepeat.ctea <- function(ctea, Alength = 10, percentA = 0.8) {
  ctea$polyA <- ctea[, grepl(paste0("A{", Alength, "}"), seq) |
                       stringr::str_count(seq, "A") / nchar(seq) > percentA]
  ctea$polyT <- ctea[, grepl(paste0("T{", Alength, "}"), seq) |
                       stringr::str_count(seq, "T") / nchar(seq) > percentA]
  ctea$normalA <- ctea[, (ori == "f" & polyT) | (ori == "r" & polyA)]
  ctea <- ctea[normalA == F, !"normalA"]  # excluding normal transcription stop site
  # TEI-polyA
  ctea[, isPolyA := rep(NA, .N)]
  ctea[ori == "f", isPolyA := polyA]
  ctea[ori == "r", isPolyA := polyT]
  ctea <- ctea[isPolyA | !grepl("^[(ACGT)]+n$", family), !c("polyA", "polyT")]
  ctea[isPolyA | !grepl("^([CT]+|[GA]+)$", seq), ]  # excluding CT or GA repeat
}

grLocator <- function(chr, start, end = start,
                      refgr, refdt = data.table(data.frame(mcols(refgr)))
                      ) {
  stopifnot(length(chr) == length(start))
  if(length(chr) == 0)
    return(refdt[0, ])
  if(missing(refgr)) {
    refgr <- end
    end <- start
  }
  isNA <- is.na(chr) | is.na(start) | is.na(end)

  ch <- chr[!isNA]
  if(grepl("chr", seqlevels(refgr)[1])) {
    ch <- pastechr(ch)
  }
  st <- pmin(start[!isNA], end[!isNA])
  en <- pmax(start[!isNA], end[!isNA])
  cgr <- GRanges(ch, IRanges(st, en), "*")
  overlap <- findOverlaps(cgr, refgr, ignore.strand = T)
  matchdata <- data.table(idx = queryHits(overlap),
                          refdt[subjectHits(overlap), ]
  )
  matchdata <- matchdata[, lapply(.SD, function(x) str_c(unique(x), collapse = ",")), by = idx]
  anno <- matchdata[match(seq_len(sum(!isNA)), idx), !"idx"]  # include non-matched position as NA
  allidx <- rep(NA_integer_, length(start))
  allidx[!isNA] <- seq_len(sum(!isNA))
  anno[allidx, ]
}

polyTElocation.ctea <- function(ctea, clip_width = 5) {
  require(GenomicRanges)
  tegr <- readRDS(polyTE_file)
  seqlevels(tegr) %<>% pastechr
  tegr$class <- NULL
  tegr$subfamily <- NULL
  polyTE <- ctea[, grLocator(chr, pos - clip_width, pos + clip_width, tegr)]
  ctea$polyTE <- polyTE$family
  ctea
}

repeatLocation <- function(chr, start, end = start) {
  # require(rtracklayer)
  # rmsk <- import(rmsk_gtf_file, "gtf")
  rmsk <- readRDS(rmsk_file)
  rmsk$family[is.na(rmsk$family)] <- rmsk$class[is.na(rmsk$family)]
  mcols(rmsk) %<>% with(DataFrame(repFamily = family, repName = subfamily))
  grLocator(chr, start, end, rmsk)
}

#' @export
repeatPositon.ctea <- function(ctea, clip_width = 5L) {
  rpt <- ctea[, repeatLocation(chr, pos - clip_width, pos + clip_width)]
  ctea[, c("posRepFamily", "posRep") := rpt]
}

#' @export
filterSimpleSite <- function(ctea,
                             threads = getOption("mc.cores", getDTthreads()),
                             Alength = 10, percentA = 0.8,
                             refAlength = 20, refPercentA = 0.8,
                             matchscore_cutoff = 0) {
  library(stringr)
  library(Biostrings)
  require(bsgenomePKG,  character.only = T)
  genome <- get(bsgenomePKG)
  
  popt <- options(mc.cores = threads)
  on.exit(options(popt))
  
  if(!exists("posRep", ctea)) {
    ctea <- repeatPosition.ctea(ctea)
  }
  
  st <- ctea[, pmax(ifelse(ori == "f", pos + 1, pos - refAlength),
                    0
  )
  ]
  en <- ctea[, pmin(ifelse(ori == "f", pos + refAlength, pos - 1),
                    seqlengths(genome)[pastechr(chr)]
  )
  ]
  refrootseq <- ctea[, getSeq(genome, pastechr(chr), st, en)] %>% as.character
  
  Asite <- stringr::str_count(refrootseq, "A") / refAlength >= refPercentA
  Tsite <- stringr::str_count(refrootseq, "T") / refAlength >= refPercentA
  
  proxseq <- ctea[, ifelse(ori == "f",
                           stringr::str_sub(seq, start = -Alength),
                           stringr::str_sub(seq, end = Alength)
  )]
  Aprop <- stringr::str_count(proxseq, "A") / Alength
  Tprop <- stringr::str_count(proxseq, "T") / Alength
  simpleSite <- (Asite & Aprop >= percentA) | (Tsite & Tprop >= percentA)
  ctea <- ctea[simpleSite == F]
  
  ctea[posRepFamily == "Simple_repeat", repeatseq := str_extract(posRep, "(?<=[(])[ATGC]+(?=[)]n)")]
  ctea[posRepFamily == "Simple_repeat", seqtomatch := str_dup(repeatseq, nchar(seq) / nchar(repeatseq))]
  ctea[posRepFamily == "Simple_repeat",
       matchscore := mcmapply(pairwiseAlignment,
                              seqtomatch,
                              seq,
                              type = "global-local") %>% mclapply(score) %>% unlist]
  ctea[is.na(matchscore) | matchscore < matchscore_cutoff, !c("repeatseq", "seqtomatch", "matchscore")]
}

filterTEsite.ctea <- function(ctea, clip_width = 5,
                              filterFamily = c("Alu", "Simple_repeat", "Low_complexity")
                     ) {
  warning("This function is deprecated. Use repeatPositon.ctea and filterSimpleSite instead.")
  rpt <- ctea[, repeatLocation(chr, pos - clip_width, pos + clip_width)]
  ctea[, c("posRepFamily", "posRep") := rpt]
  ctea <- ctea[strsplit(posRepFamily, ",") %>%
         sapply(function(x) ! any(x %in% filterFamily)), ]
  ctea[!grepl("SVA", family) | !grepl("SVA", posRep)]
}

filterPolyTEsite.ctea <- function(ctea, clip_width = 5) {
  warning("This function is deprecated. Use polyTElocation.ctea instead.")
  require(GenomicRanges)
  tegr <- readRDS(polyTE_file) %>%
  {resize(., width(.) + clip_width, fix = "start")} %>%
  {resize(., width(.) + clip_width, fix = "end")}
  cgr <- ctea[, GRanges(pastechr(chr),
                        IRanges(pos, width = 1),
                        strand = rep("*", .N))
              ]
  overlapTE <- findOverlaps(cgr, tegr, ignore.strand = T)
  ctea[, inPolyTE := .I %in% queryHits(overlapTE)]
  ctea[inPolyTE == F, ]
}

filterIn.ctea <- function(ctea, xtea, accept_width = 3) {
  require(GenomicRanges)
  cgr <- ctea[, GRanges(chr,
                        IRanges(pos, width = 1),
                        strand = c(r = "+", f = "-")[ori]
  )
  ]
  xgr <- xtea[, GRanges(chr,
                        IRanges(pos, width = 1),
                        strand = c(r = "+", f = "-", "*" = "*")[ori]
  ) %>% flank(accept_width, start = F, both = T)
  ]
  overlaps <- findOverlaps(cgr, xgr)
  ctea[, inxtea := .I %in% queryHits(overlaps)]
  ctea[inxtea == T, ]
}

filterNoClip.ctea <- function(ctea, similarity_cutoff = 75, threads = getOption("mc.cores", detectCores())) {
  require(bsgenomePKG,  character.only = T)
  genome <- get(bsgenomePKG)
  ctea[, start := pos]
  ctea[, end := pos]
  if(nrow(ctea) == 0) {
    msg("No insertion found.")
    ctea[, similar := numeric(0)]
    return(ctea)
  }
  ctea[ori == "f", start := pos - nchar(seq) + 1L]
  ctea[ori == "r", end := pos + nchar(seq) - 1L]
  overstart <- ctea[, pmax(0, 1 - start)]
  overend <- ctea[, pmax(0, end - seqlengths(genome)[pastechr(chr)])]
  ctea[, start := start + overstart]
  ctea[, end := end - overend]
  refseq <- ctea[, getSeq(genome, pastechr(chr), start, end)]
  refseq <- paste0(stringr::str_dup("N", overstart), refseq, stringr::str_dup("N", overend))
  n <- nrow(ctea)
  msg("Comparing clipped sequences with reference sequences...")
  similar <- mclapply(
    seq_len(n),
    mc.cores = threads,
    mc.preschedule = T,
    function(i) {
      cat(sprintf("\r%0.2f%%              ", i/n*100))
      rs <- refseq[i]
      cs <- ctea[i, seq]
      pid(pairwiseAlignment(rs, cs))
    }
  )

  msg("Comparing done.")
  ctea[, similar := unlist(similar)]
  ctea[similar <= similarity_cutoff, ]
}

checkStrandedness <- function(bamfile, edbpkg = edbPKG) {
  msg("Checking strandedness of bamfile...")

  library(Rsamtools)
  require(edbpkg, character.only = T)
  edb <- get(edbpkg)

  housekeeping <- c("RRN18S", "GAPDH", "PGK1", "PPIA", "RPLP0",
                    "ARBP", "B2M", "YWHAZ", "SDHA", "TFRC",
                    "GUSB", "HMBS", "HPRT1", "TBP"
  )

  he <- exons(edb, filter = GeneNameFilter(housekeeping), column = "tx_biotype")
  he %<>% subset(tx_biotype %in% "protein_coding")

  nohe <- exonsByOverlaps(edb, he, columns = c("gene_name")) %>%
    subset(!gene_name %in% housekeeping)
  ovlgene <- subsetByOverlaps(he, nohe, ignore.strand = T)$gene_name %>% unique
  he %<>% subset(!gene_name %in% ovlgene)

  bamheader <- scanBamHeader(bamfile)
  seqlevels(he, pruning.mode = "coarse") <- intersect(seqlevels(he), names(bamheader[[1]]$targets))

  readcnt <- countBam(bamfile, param = ScanBamParam(which = he)) %>% data.table
  hemc <- readcnt[, records %between% quantile(records, c(.45, .55))] %>% he[.]
  sam <- scanBam(bamfile,
                 param = ScanBamParam(which = hemc,
                                      what = c("flag", "strand"),
                                      simpleCigar = T,
                                      mapqFilter = 30L
                 )
  )
  isfirst <- sapply(sam, `[[`, "flag") %>% unlist %>% bamFlagTest("isFirstMateRead")
  rstrand <- sapply(sam, `[[`, "strand") %>% unlist
  gstrand <- rep(strand(hemc), sapply(sam, function(x) length(x[[1]]))) %>% as.factor
  tbl <- table(isfirst, rstrand == gstrand)
  
  if(dim(tbl)[1] < 2 | dim(tbl)[2] < 2) {
    return("non-stranded")
  } else {
    stat <- fisher.test(tbl)
    stranded <- stat$p.value < .05 & abs(log2(stat$estimate)) > 2
    if(!stranded) {
      "non-stranded"
    } else {
      if(mean(isfirst[rstrand == gstrand]) > 0.5) {
        "first-sense"
      } else {
        "first-antisense"
      }
    }
  }

}

subsampleBam <- function(bamfile,
                         gr,
                         numReads = countBam(
                           bamfile,
                           param = ScanBamParam(which = gr, mapqFilter = mapqFilter)
                         )$records,
                         what = c("qname", "seq", "qual", "mapq", "flag", "mpos"),
                         tag = "NM",  # character(0)
                         mapqFilter = 1L,
                         maxReads = 1e5) {

  gr %<>% reduce
  cbfile <- tempfile()
  on.exit(unlink(cbfile))
  if(numReads <= maxReads) {
    sam <- readGAlignments(
      bamfile,
      param = ScanBamParam(
        which = gr,
        what = what,
        tag = tag,
        mapqFilter = mapqFilter
      )
    )
    return(sam)
  }
  region <- paste0(seqnames(gr), ":", start(gr), "-", end(gr), collapse = " ")
  msg("Subsampling ", region, " from ", numReads, " reads to ", maxReads)
  awk <- "'BEGIN {srand()} /^@/ {print} !/^@/ {if(rand() * n-- < p) {p--; print; if(p==0) exit}}'"
  cmd <- paste("samtools view -h -q", mapqFilter, bamfile, region, "|",
               "awk", awk, sprintf("n=%d p=%d", numReads, maxReads), "|",
               "samtools view -bh -", ">",
               cbfile
               )
  # cmd <- paste("samtools view -h -q", mapqFilter, bamfile, region, "|",
  #              "Rscript --vanilla", subsample_script, numreads, maxReads, "|",
  #              "samtools view -b -", ">",
  #              cbfile)
  # writeLines(cmd)
  exitcode <- system(cmd)
  if(exitcode != 0) {
    stop("Error while subsampling bam file: ", region)
  }
  subsamplecnt <- countBam(cbfile)$records
  if(subsamplecnt != maxReads) {
    stop("Error: subsampling incomplete in ", region)
  }
  readGAlignments(cbfile,
                  param = ScanBamParam(
                    what = what,
                    tag = tag,
                    mapqFilter = mapqFilter
                  ))
}


getClippedReads <- function(bamfile, chr, pos, ori = c("f", "r"),
                            searchWidth = 10L,
                            mapqFilter = 1L,
                            subsample = F,
                            numReads = NULL,
                            maxReads = 1e5) {
  require(GenomicAlignments)
  gr <- GRanges(chr,
                IRanges(min(pos) - searchWidth, max(pos) + searchWidth),
                strand="*"
  )

  sam <- if(subsample && numReads > maxReads) {
    subsampleBam(bamfile,
                 gr,
                 numReads = numReads,
                 what = c("qname", "seq", "qual", "mapq", "flag", "mpos"),
                 tag = "NM",
                 mapqFilter = mapqFilter,
                 maxReads = maxReads
    )
  } else {
    readGAlignments(
      bamfile,
      param = ScanBamParam(
        which = gr,
        what = c("qname", "seq", "qual", "mapq", "flag", "mpos"),
        tag = "NM",
        mapqFilter = mapqFilter
      )
    )
  }

  assignZero <- function(mcols) {
    mcols$sseq <- character(0)
    mcols$squal <- numeric(0)
    mcols$overhang <- integer(0)
    mcols$gap <- integer(0)
    mcols
  }

  if(length(sam) == 0) {
    mcols(sam) <- assignZero(mcols(sam))
    return(sam)
  }

  if(ori == "f") {
    # deduplication
    sam <- sam[start(sam) > pos - qwidth(sam)]
    secondary <- bamFlagTest(mcols(sam)$flag, "isSecondaryAlignment")
    split <- grepl("^\\d+S", cigar(sam))
    othersplit <- endsWith(cigar(sam), "S")
    sam <- sam[order(secondary, split, othersplit)]
    isFirst <- bamFlagTest(mcols(sam)$flag, "isFirstMateRead")
    sam <- sam[!duplicated(paste(mcols(sam)$qname, isFirst)) & grepl("^\\d+S", cigar(sam))]
    # mcols(sam)$qname  <- NULL

    mcols(sam)$slen <- sub("S.*", "", cigar(sam)) %>% as.integer
    mcols(sam)$shift <- start(sam) - pos - 1
    if(length(sam) == 0) {
      mcols(sam) <- assignZero(mcols(sam))
      return(sam)
    }
    sam <- sam[abs(mcols(sam)$shift) < searchWidth &
                 mcols(sam)$slen - mcols(sam)$shift > 0 &
                 mcols(sam)$slen - mcols(sam)$shift <= width(mcols(sam)$seq)]
    mcols(sam)$sseq <- subseq(mcols(sam)$seq, 1, mcols(sam)$slen - mcols(sam)$shift)
    mcols(sam)$squal <- subseq(mcols(sam)$qual, 1, mcols(sam)$slen - mcols(sam)$shift) %>%
      as("IntegerList") %>%
      mean
    mcols(sam)$qual <- NULL
    mcols(sam)$overhang <- stringr::str_extract(cigar(sam), "(?<=S)[0-9]+(?=M[0-9]+N)") %>%
      as.integer
    mcols(sam)$gap <- stringr::str_extract(cigar(sam), "[0-9]+(?=N)") %>%
      as.integer

  } else if(ori == "r") {
    # deduplication
    sam <- sam[end(sam) < pos + qwidth(sam)]
    secondary <- bamFlagTest(mcols(sam)$flag, "isSecondaryAlignment")
    split <- endsWith(cigar(sam), "S")
    othersplit <- grepl("^\\d+S", cigar(sam))
    sam <- sam[order(secondary, split, othersplit)]
    isFirst <- bamFlagTest(mcols(sam)$flag, "isFirstMateRead")
    sam <- sam[!duplicated(paste(mcols(sam)$qname, isFirst)) & endsWith(cigar(sam), "S")]
    # mcols(sam)$qname  <- NULL

    mcols(sam)$slen <-
      explodeCigarOpLengths(cigar(sam), ops="S") %>%
      sapply(function(x) x[length(x)])
    mcols(sam)$shift <- end(sam) - pos + 1
    if(length(sam) == 0) {
      mcols(sam) <- assignZero(mcols(sam))
      return(sam)
    }
    sam <- sam[abs(mcols(sam)$shift) < searchWidth &
                 mcols(sam)$slen + mcols(sam)$shift > 0 &
                 mcols(sam)$slen + mcols(sam)$shift <= width(mcols(sam)$seq)]
    mcols(sam)$sseq <- subseq(mcols(sam)$seq, -(mcols(sam)$slen + mcols(sam)$shift))
    mcols(sam)$squal <- subseq(mcols(sam)$qual, -(mcols(sam)$slen + mcols(sam)$shift)) %>%
      as("IntegerList") %>%
      mean
    mcols(sam)$qual <- NULL
    mcols(sam)$overhang <- stringr::str_extract(cigar(sam), "(?<=N)[0-9]+(?=M[0-9]+S)") %>%
      as.integer
    mcols(sam)$gap <- stringr::str_extract(cigar(sam), "[0-9]+(?=N[^N]+$)") %>%
      as.integer
  }

  sam
}

#' @param gal \code{GAlignments} object
isposclipped <- function(gal, ori, pos, searchWidth = 10) {
  clipped <- grepl(ifelse(ori == "f", "^[0-9]+S", "S$"),
                   cigar(gal)
  )
  cigarpos <- if(ori == "f") {
    start(gal) - 1
  } else {
    # cigarrle <- cigarToRleList(cigar(gal))
    # start(gal) +
    #   sum(cigarrle %in% c("M", "N", "D"))
    end(gal) + 1
  }
  clipped & cigarpos %between% c(pos - searchWidth, pos + searchWidth)
}

getClippedPairs <- function(bamfile, chr, pos, ori = c("f", "r"),
                            firstsam = NULL,
                            what = c("qname", "flag"),
                            tag = character(0),
                            noOverClip = T,
                            searchWidth = 10L,
                            mapqFilter = 1L,
                            subsample = F,
                            numReads = NULL,
                            maxReads = 1e6) {
  require(GenomicAlignments)

  if(missing(firstsam)) {
    gr <- GRanges(chr,
                  IRanges(min(pos) - searchWidth, max(pos) + searchWidth),
                  strand="*"
    )

    if(subsample && missing(numReads)) {
      numReads <- countBam(bamfile,
                           param = ScanBamParam(which = gr, mapqFilter = mapqFilter)
      )$records
    }
    sam <- if(subsample && numReads > maxReads) {
      subsampleBam(bamfile, gr,
                   what = c(what, "mpos"),
                   tag = tag,
                   mapqFilter = mapqFilter,
                   numReads = numReads,
                   maxReads = maxReads
      )
    } else {
      readGAlignments(bamfile,
                      param = ScanBamParam(
                        which = gr,
                        what = c(what, "mpos"),
                        tag = tag,
                        mapqFilter = mapqFilter
                      )
      )
    }

    sam %<>% .[isposclipped(sam, ori, pos, searchWidth = searchWidth)]
  } else {
    sam <- firstsam
  }

  sam %<>% .[order(njunc(.), decreasing = T)]
  mgr <- with(subset(sam, bamFlagTest(flag, "isProperPair")),
       GRanges(seqnames, IRanges(mpos, mpos), invertStrand(strand))
  )
  if(noOverClip) {
    mgr <- if(ori == "f") {
      mgr[start(mgr) >= pos]
    } else {
      mgr[end(mgr) <= pos]
    }
  }
  mgr %<>% resize(100) %>% intersect(range(mgr))

  if(length(mgr) == 0) {
    return(list(
      pairs = GAlignmentPairs(sam[0], sam[0]),
      nopair = sam
    ))
  }

  mrecords <- countBam(
    bamfile,
    param = ScanBamParam(which = mgr,
                         mapqFilter = mapqFilter)
  )$records
  mate <- if(subsample && mrecords > maxReads) {
    subsampleBam(bamfile, mgr,
                 numReads = mrecords,
                 what = what,
                 tag = tag,
                 mapqFilter = mapqFilter,
                 maxReads = maxReads
    )
  } else {
    readGAlignments(bamfile,
                    param = ScanBamParam(
                      which = mgr,
                      what = what,
                      tag = tag,
                      mapqFilter = mapqFilter
                    )
    )
  }
  mate %<>% subset(qname %in% mcols(sam)$qname)
  mate %<>% .[order(njunc(.), decreasing = T)]

  sam %<>% .[!duplicated(mcols(.)$qname)]
  sam1 <- subset(sam, bamFlagTest(flag, "isFirstMateRead"))
  sam2 <- subset(sam, !bamFlagTest(flag, "isFirstMateRead"))
  mate1 <- subset(mate, bamFlagTest(flag, "isFirstMateRead"))
  mate2 <- subset(mate, !bamFlagTest(flag, "isFirstMateRead"))
  m1 <- match(mcols(sam1)$qname, mcols(mate2)$qname)
  m2 <- match(mcols(sam2)$qname, mcols(mate1)$qname)
  pairs <- GAlignmentPairs(c(sam1[!is.na(m1)], sam2[!is.na(m2)]),
                           c(mate2[na.omit(m1)], mate1[na.omit(m2)]))
  nopair <- c(sam1[is.na(m1)], sam2[is.na(m2)])
  list(pairs = pairs, nopair = nopair)
}

compareClippedSeq <- function(sseq, seq, ori,
                              shift_range = 0, shift_dir = "both",
                              min_length = 1) {
  require(Biostrings)

  if(ori == "f") {
    subseqprox <- function(s, n) subseq(s, -n)
    shiftseq <- function(s) subseq(s, 1, -2)
  } else if(ori == "r") {
    subseqprox <- function(s, n) subseq(s, 1, n)
    shiftseq <- function(s) subseq(s, 2)
  }

  longseq <- nchar(sseq) > nchar(seq)
  sseq[longseq] <- subseqprox(sseq[longseq], nchar(seq))
  sseq[nchar(sseq) < min_length] <- "N"
  gseq <- if(ori == "f") {
    paste0(stringr::str_dup("N", nchar(seq) - nchar(sseq)), sseq)
  } else {
    paste0(sseq, stringr::str_dup("N", nchar(seq) - nchar(sseq)))
  }
  
  seqm <- strsplit(gseq, "") %>% do.call(rbind, .)
  seqm[seqm == "N"] <- NA_character_
  differ <- sweep(seqm, 2, strsplit(seq, "")[[1]], "!=") %>% rowSums(na.rm = T) %>% {. / nchar(sseq)}
  if(nchar(seq) < 2) shift_range <- 0
  if(shift_range <= 0) {
    differ
  } else if(shift_dir == "both") {
    pmin(
      differ,
      compareClippedSeq(sseq, shiftseq(seq), ori,
                        shift_range = shift_range - 1,
                        shift_dir = "dist",
                        min_length = min_length + 1),
      compareClippedSeq(shiftseq(sseq), seq, ori,
                        shift_range = shift_range - 1,
                        shift_dir = "prox",
                        min_length = min_length + 1),
      na.rm = T
    )
  } else {
    if(shift_dir == "dist") {
      seq <- shiftseq(seq)
    } else if(shift_dir == "prox") {
      sseq <- shiftseq(sseq)
    } else {
      stop('The value of "shift_dir" must be either "prox" or "dist"')
    }
    pmin(
      differ,
      compareClippedSeq(sseq, seq, ori,
                        shift_range = shift_range - 1,
                        shift_dir = shift_dir,
                        min_length = min_length + 1)
    )
  }
}


consensusClip <- function(bamfile, chr, pos, ori, searchWidth = 10) {
  require(Biostrings)
  stopifnot(length(pos) == 1)
  if(is.na(pos)) return(NA)
  sam <- getClippedReads(bamfile, chr, pos, ori, searchWidth)
  sseq <- mcols(sam)$sseq[mcols(sam)$squal > 20]
  if(length(sseq) == 0) return("")
  if(ori == "f") sseq <- reverse(sseq)
  consensus <- consensusString(sseq, ambiguityMap = "N", threshold = 0.5)
  if(ori == "f") consensus <- reverse(consensus)
  if(nchar(consensus) == 0) {
    NA
  } else {
    as.character(consensus)
  }
}

TEalignScore <- function(seq, TEfamily) {
  if(length(seq) == 0) return(numeric(0))
  stopifnot(length(TEfamily) == 1)
  library(Biostrings)
  fa <- readDNAStringSet(refTEfa)
  names(fa) %<>% sub("\\s.*", "", .)
  if(!TEfamily %in% names(fa)) return(NA)
  fa %<>% .[[TEfamily]]
  fscore <- pairwiseAlignment(seq, fa, type = "global-local") %>%
    score
  rscore <- DNAStringSet(seq) %>%
    reverseComplement %>%
    pairwiseAlignment(fa, type = "global-local") %>%
    score
  if(sum(fscore) > sum(rscore)) {
    fscore
  } else {
    rscore
  }
}

cntFilter.ctea <- function(ctea,
                           trueCntCutoff = 2,
                           uniqueCntCutoff = 2,
                           falseCntRatioCutoff = 10,
                           baseQualCutoff = 25,
                           wrongPosPropCutoff = 0.2,
                           bothClipPropCutoff = 0.4,
                           secondaryCutoff = 0.99,
                           lowMapQualPropCutoff = 0.99,
                           TEscoreCutoff = 30,
                           nonspecificTEcutoff = 0,
                           hardFilter_cutoff = 50) {
  if(!exists("Filter", ctea)) {
    ctea[, Filter := ""]
  } else {
    ctea[Filter == "PASS", Filter := ""]
  }
  ctea[trueCnt <= trueCntCutoff | uniqueCnt <= uniqueCntCutoff,
       Filter := paste(Filter, "lowCnt", sep=";")]
  ctea[falseCnt / trueCnt >= falseCntRatioCutoff, Filter := paste(Filter, "noisy", sep=";")]
  ctea[baseQual <= baseQualCutoff, Filter := paste(Filter, "lowPhred", sep=";")]
  ctea[anyWrongPos / matchCnt >= wrongPosPropCutoff, Filter := paste(Filter, "badMap", sep=";")]
  ctea[bothClip / matchCnt > bothClipPropCutoff, Filter := paste(Filter, "pseudoMap", sep=";")]
  ctea[secondary >= secondaryCutoff, Filter := paste(Filter, "secondary", sep=";")]
  ctea[lowMapQual / trueCnt >= lowMapQualPropCutoff, Filter := paste(Filter, "lowMapQual", sep=";")]
  if(exists("TEscore", ctea)) {
    ctea[TEscore <= TEscoreCutoff, Filter := paste(Filter, "lowTEscore", sep = ";")]
  }
  if(exists("nonspecificTE", ctea)) {
    ctea[nonspecificTE >= nonspecificTEcutoff, Filter := paste(Filter, "nonspecificTE", sep = ";")]
  }
  if(exists("hardDist", ctea)) {
    ctea[hardDist < hardFilter_cutoff, Filter := paste(Filter, "indel", sep=";")]
  }
  if(exists("hardSpl", ctea) & exists("hardTE", ctea)) {
    ctea[!is.na(hardSpl) & !is.na(hardTE), Filter := paste(Filter, "notGapped", sep=";")]
  }
  ctea[, Filter := sub("^;", "", Filter)]
  ctea[Filter == "", Filter := "PASS"]
  ctea[, Filter := sapply(strsplit(Filter, ";"),
                          function(x) paste(unique(x), collapse = ";"))]
  ctea
}

countClippedReads.ctea <- function(ctea,
                                   bamfile,
                                   edbpkg = edbPKG,
                                   strandedness = c("non-stranded", "first-sense", "first-antisense"),
                                   searchWidth = 10L,
                                   mapqFilter = 1L,
                                   shift_range = 0,
                                   mismatch_cutoff = 0.1,
                                   cliplength_cutoff = 4,
                                   maxReads = 1e5,
                                   fusiontype_cutoff = 3,
                                   threads = getOption("mc.cores", detectCores())) {
  library(BiocParallel)
  require(GenomicAlignments)
  require(data.table)

  NAresult <- list(
    depth = NA_real_,
    matchCnt = NA_integer_,
    trueCnt = NA_integer_,
    uniqueCnt = NA_integer_,
    bothClip = NA_integer_,
    polyAcnt = NA_integer_,
    falseCnt = NA_integer_,
    discCnt = NA_integer_,
    baseQual = NA_real_,
    lowMapQual = NA_integer_,
    clustered = NA_real_,
    anyOverClip = NA_integer_,
    mateDist = NA_integer_,
    anyWrongPos = NA_integer_,
    overhang = NA_integer_,
    gap = NA_integer_,
    secondary = NA_real_,
    editDistance = NA_real_,
    nonspecificTE = NA_real_,
    r1pstrand = NA_real_,
    fusion_tx_id = NA_character_,
    tx_support_exon = NA_integer_,
    tx_support_intron = NA_integer_,
    numgap = NA_integer_
  )

  n <- nrow(ctea)
  if(n == 0) {
    cntdt <- as.data.table(NAresult)[0]
    return( data.table(ctea, cntdt) )
  }

  if(missing(strandedness)) {
    strandedness <- checkStrandedness(bamfile, edbpkg)
  }
  if(paste0("package:", edbpkg) %in% search()) {
    detach(paste0("package:", edbpkg), unload = T, character.only = T)
  }

  msg("Counting clipped reads...")
  countThis <- function(i, bamfile, chr, pos, ori, seq, family) {
    cat(sprintf("\r%0.2f%%              ", i/n*100))
    if(anyNA(c(chr, pos, ori, seq))) {
      return(NAresult)
    }
    records <- countBam(
      bamfile,
      param = ScanBamParam(which = GRanges(chr, IRanges(pos - searchWidth, pos + searchWidth), "*"),
                           mapqFilter = mapqFilter)
    )$records
    sam <- getClippedReads(bamfile,
                           chr, pos, ori,
                           subsample = records > maxReads,
                           searchWidth = searchWidth,
                           mapqFilter = mapqFilter,
                           numReads = records,
                           maxReads = maxReads
    )

    if(length(sam) == 0) {
      res <- NAresult
      res$depth <- records / (2 * searchWidth + 1)
      cntidx <- grep("Cnt", names(res))
      res[cntidx] <- 0
      return(res)
    }

    meta <- mcols(sam)
    differ <- compareClippedSeq(meta$sseq, seq, ori, shift_range = shift_range)
    isMatch <- differ < mismatch_cutoff
    bothClip <- grepl("^[0-9]+S.*S$", cigar(sam))
    shortClip <- nchar(meta$sseq) < cliplength_cutoff
    isPolyA <- grepl("^(A+|T+)$", meta$sseq)
    isMateSide <- as.factor(strand(sam)) == c(f = "-", r = "+")[ori]
    isProperPair <- bamFlagTest(meta$flag, "isProperPair")
    mateUnmapped <- bamFlagTest(meta$flag, "hasUnmappedMate")
    isSecondary <- bamFlagTest(meta$flag, "isSecondaryAlignment")
    mateOverlap <- meta$mpos %between% list(start(sam) - 3, end(sam) + 3)
    isOverClip <- isProperPair & !mateOverlap & isMateSide & isMatch & !isPolyA
    possibleDup <- duplicated(paste(nchar(meta$sseq), meta$mpos)) | duplicated(meta$sseq)
    isTEread <- TEalignScore(meta$seq, family)
    r1strand <- ifelse(bamFlagTest(meta$flag, "isFirstMateRead"), strand(sam), invertStrand(strand(sam)))

    cnt <- list(
      depth = records / (2 * searchWidth + 1),
      matchCnt = sum(isMatch),
      trueCnt = sum(isMatch & !shortClip & !isPolyA),
      uniqueCnt = sum(isMatch & !possibleDup & !shortClip),
      bothClip = sum(isMatch & bothClip),
      polyAcnt = sum(isPolyA),
      falseCnt = sum(!isMatch & !shortClip),
      discCnt = sum(isMatch & !isProperPair & !mateUnmapped & isMateSide),
      baseQual = median(meta$squal[isMatch & !isPolyA]),
      lowMapQual = sum(meta$mapq[isMatch & !shortClip] < mapqFilter + 10),
      clustered = sd(nchar(meta$sseq[isMatch])),
      anyOverClip = sum(isOverClip),
      mateDist = suppressWarnings(min(abs(pos - meta$mpos)[isOverClip])),
      anyWrongPos = sum(!isProperPair & !mateUnmapped & !isMateSide & isMatch),
      overhang = suppressWarnings(median(meta$overhang[isMatch & !shortClip])),
      gap = suppressWarnings(getmode(meta$gap[isMatch & !shortClip])),
      secondary = mean(isSecondary[isMatch]),
      editDistance = mean(meta$NM[isMatch]),
      nonspecificTE = mean(isTEread[isMatch]),
      # fusion_tx_id = fusionTx$tx_id,
      # tx_support_exon = fusionTx$tx_support_exon,
      # tx_support_intron = fusionTx$tx_support_intron,
      # numgap = fusionTx$numgap,
      r1pstrand = mean(r1strand[isMatch] == "+")
    )

    if(sum(isMatch) >= fusiontype_cutoff) {

      mcols(sam)$seq <- NULL
      mcols(sam)$sseq <- NULL
      sampair <- getClippedPairs(bamfile, chr, pos, ori, sam[isMatch],
                                 subsample = T, maxReads = maxReads) %>% {
        c(first(.$pairs), second(.$pairs), .$nopair)
      }
      fusionTx <- fusionTxMatch(sampair, edbpkg, strandedness)
    } else {
      fusionTx <- data.table(tx_id = NA_character_,
                             tx_support_exon = NA_integer_,
                             tx_support_intron = NA_integer_,
                             numgap = NA_integer_
      )
    }
    setnames(fusionTx, "tx_id", "fusion_tx_id")

    c(cnt, as.list(fusionTx))

  }
  lcnt <- bptry(
    bpmapply(
      FUN = countThis,
      seq_len(n),
      ctea$chr, ctea$pos, ctea$ori, ctea$seq, ctea$family,
      MoreArgs = list(bamfile = bamfile),
      SIMPLIFY = F,
      BPPARAM = MulticoreParam(workers = threads, stop.on.error = T, tasks = n)
    ),
    bplist_error = function(e) {
      message(e, "\n")
      save(ctea, bamfile,
           searchWidth, mapqFilter, shift_range, mismatch_cutoff, cliplength_cutoff, maxReads, threads,
           e,
           file = "countClippedReadsErr.RData")
      e
      # quit("no", 1)
    }
  )

  msg("Counting done.")
  cntdt <- tryCatch(
    rbindlist(lcnt),
    error = function(e) {
      message(e, "\n")
      save(ctea, bamfile,
           searchWidth, mapqFilter, shift_range, mismatch_cutoff, cliplength_cutoff, maxReads, threads,
           lcnt,
           e,
           file = "countClippedReadsRbindErr.RData")
      e
      # quit("no", 1)
    }
  )
  stopifnot(ctea[, .N] == cntdt[, .N])

  cntdt[, strand := "*"]
  if(strandedness == "first-sense") {
    cntdt[r1pstrand > .5, strand := "+"]
    cntdt[r1pstrand < .5, strand := "-"]
  }
  if(strandedness == "first-antisense") {
    cntdt[r1pstrand > .5, strand := "-"]
    cntdt[r1pstrand < .5, strand := "+"]
  }

  data.table(ctea, cntdt)
}

countClippedReads.default <- function(chr, pos, ori, seq,
                                      bamfile,
                                      searchWidth = 10,
                                      shift_range = 0,
                                      mismatch_cutoff = 0.1,
                                      cliplength_cutoff = 4,
                                      threads = getOption("mc.cores", detectCores())) {
  require(data.table)
  ctea <- data.table(chr = chr, pos = pos, ori = ori, seq = seq)
  countClippedReads.ctea(
    ctea,
    bamfile,
    searchWidth = searchWidth,
    shift_range = shift_range,
    mismatch_cutoff = mismatch_cutoff,
    cliplength_cutoff = cliplength_cutoff,
    threads = threads
  )
}

countBothClippedReads <- function(chr,
                                  fpos, fseq,
                                  rpos, rseq,
                                  bamfile,
                                  searchWidth = 10,
                                  shift_range = 0,
                                  mismatch_cutoff = 0.1,
                                  cliplength_cutoff = 4,
                                  threads = getOption("mc.cores", detectCores())) {
  require(data.table)
  fcnt <- countClippedReads.default(
    chr, fpos, "f", fseq,
    bamfile,
    searchWidth = searchWidth,
    shift_range = shift_range,
    mismatch_cutoff = mismatch_cutoff,
    cliplength_cutoff = cliplength_cutoff,
    threads = threads
  )
  fcnt[, ori := NULL]
  names(fcnt)[-1] <- paste0("f", names(fcnt)[-1])
  rcnt <- countClippedReads.default(
    chr, rpos, "r", rseq,
    bamfile,
    searchWidth = searchWidth,
    mismatch_cutoff = mismatch_cutoff,
    cliplength_cutoff = cliplength_cutoff,
    threads = threads
  )
  rcnt[, ori := NULL]
  names(rcnt)[-1] <- paste0("r", names(rcnt)[-1])
  data.table(fcnt, rcnt[, !"chr"])
}

ungapPos.rtea <- function(rtea, overhang_cutoff = 5L) {
  ungappos <- rtea$pos
  fover <- rtea[, which(ori == "f" & overhang <= overhang_cutoff)]
  ungappos[fover] <- rtea[fover, pos + gap]
  rover <- rtea[, which(ori == "r" & overhang <= overhang_cutoff)]
  ungappos[rover] <- rtea[rover, pos - gap]
  ungappos
}

calculateDepth.default <- function(bamfile, chr, pos, ori, width = 10L) {
  region <- if(ori == "f") {
    sprintf("%s:%d-%d", chr, pos, pos + width)
  } else {
    sprintf("%s:%d-%d", chr, pos - width - 1, pos - 1)
  }
  samtoolsCmd <- paste("samtools", "depth",
                       "-Q", 1,
                       "-g", 4095,
                       "-a",
                       bamfile,
                       "-r", region)
  posDepth <- fread(cmd = samtoolsCmd,
                    col.names = c("chr", "pos", "depth"))
  max(posDepth$depth)
}

calculateDepth.rtea <- function(rtea, bamfile, width = 10L, threads = getOption("mc.cores", detectCores())) {
  samtoolsVer <- system("samtools --version", intern = T)
  majorVer <- sub(".* ", "", samtoolsVer[1]) %>% sub("[.].*", "", .) %>% as.integer
  minorVer <- sub(".* ", "", samtoolsVer[1]) %>% sub(".*[.]", "", .) %>% as.integer
  if(majorVer < 1 | minorVer < 10) {
    stop("samtools version 1.10 or higher is required.")
  }
  
  ungappos <- ungapPos.rtea(rtea, 5L)
  depths <- rtea[, mcmapply(calculateDepth.default, 
                            chr, ungappos, ori,
                            MoreArgs = list(bamfile = bamfile,
                                            width = width),
                            mc.cores = threads)]
  rtea$depth <- depths
  rtea
}


unique.rtea <- function(rtea, ...) {
  dup <- duplicate.rtea(rtea, ...)
  rtea[is.na(dup) | dup == seq_along(dup)]
}

duplicate.rtea <- function(rtea,
                           selectbest = T,
                           verbose = getOption("verbose", F),
                           overhang_cutoff = 10L,
                           maxgap = 5L,
                           seqmatch_shift_range = maxgap,
                           seqmatch_cutoff = 0.7,
                           threads = 1) {
  require(bsgenomePKG,  character.only = T)
  genome <- get(bsgenomePKG)
  library(GenomicRanges)
  library(Biostrings)
  library(BiocParallel)
  opt <- options(mc.cores = threads)
  prevthreads <- setDTthreads(threads)
  on.exit({
    options(opt)
    setDTthreads(prevthreads)
  })
  
  # finding positional overlap
  if(verbose) msg("Finding positional overlaps...")
  ungappos <- ungapPos.rtea(rtea, overhang_cutoff = overhang_cutoff)
  overhangint <- rtea[, floor(overhang)]
  posmatch <- function(locdt, ...) {
    uniqloc <- locdt[, .(idx = list(.I)), by = .(chr, pos, ori)]
    rgr <- uniqloc[, GRanges(chr, IRanges(pos, pos), ifelse(ori == "f", "+", "-"))]
    ovl <- findOverlaps(rgr, drop.self = T, ...)
    uniqloc[, ovlidx := mclapply(ovl, function(x) unlist(idx[x]))]
    
    filterovl <- function(idx1, idx2, verbose = getOption("verbose", F)) {
      select <- idx1 < idx2
      idx1 %<>% .[select]
      idx2 %<>% .[select]
      rm(select)
      gc()
      data.table(idx1, idx2, posgap = locdt[idx2, pos] - locdt[idx1, pos])
    }
    idx1 <- uniqloc[, rep(unlist(idx), rep(elementNROWS(ovlidx), elementNROWS(idx)))]
    idx2 <- uniqloc[, rep(ovlidx, elementNROWS(idx)) %>% unlist]
    if(is.null(idx2)) idx2 <- integer(0)
    gapovlap <- filterovl(idx1, idx2)
    
    posovlap <- filterovl(
      uniqloc[elementNROWS(idx) > 1L, rep(unlist(idx), rep(elementNROWS(idx), elementNROWS(idx)))],
      uniqloc[elementNROWS(idx) > 1L, rep(idx, elementNROWS(idx)) %>% unlist]
    )
    
    rbind(gapovlap, posovlap)
  }
  ovlap <- posmatch(rtea[, .(chr, pos, ori, seq)], maxgap = maxgap)
  ugovlap <- posmatch(data.table(rtea[, .(chr, ori, seq)], pos = ungappos), maxgap = maxgap)
  ovlap %<>% funion(ugovlap)
  
  # get base sequence
  if(verbose) msg("Getting extended sequence ...")
  baseseqdt <- data.table(rtea[, .(chr, pos, ori)], overhangint, ungappos) %>%
    .[unique(ovlap[, c(idx1, idx2)]), ] %>%
    unique
  baseseqdt[ori == "f", start := pos + 1]
  baseseqdt[ori == "f", end := pos + maxgap + 1]
  baseseqdt[ori == "r", start := pos - maxgap - 1]
  baseseqdt[ori == "r", end := pos - 1]
  baseseqdt[ori == "f" & overhangint <= maxgap + 1, end := pos + overhangint]
  baseseqdt[ori == "f" & overhangint <= maxgap + 1, start2 := ungappos + overhangint + 1]
  baseseqdt[ori == "f" & overhangint <= maxgap + 1, end2 := ungappos + maxgap + 1]
  baseseqdt[ori == "r" & overhangint <= maxgap + 1, start := pos - overhangint]
  baseseqdt[ori == "r" & overhangint <= maxgap + 1, start2 := ungappos - maxgap - 1]
  baseseqdt[ori == "r" & overhangint <= maxgap + 1, end2 := ungappos - overhangint - 1]
  baseseqdt$seq <- bpvec(
    seq_len(nrow(baseseqdt)),
    function(i) baseseqdt[i, getSeq(genome, pastechr(chr), start, end)],
    BPPARAM = MulticoreParam(workers = threads)
  ) %>% as.character
  baseseqdt[overhangint <= maxgap + 1,
            seq2 := bpvec(seq_along(start2),
                          function(i) getSeq(genome, pastechr(chr[i]), start2[i], end2[i]),
                          BPPARAM = MulticoreParam(workers = threads)
            ) %>% as.character]
  baseseqdt[ori == "f" & overhangint <= maxgap + 1, seq := paste0(seq, seq2)]
  baseseqdt[ori == "r" & overhangint <= maxgap + 1, seq := paste0(seq2, seq)]
  baseseqdt[, c("start", "end", "start2", "end2", "seq2") := NULL]
  
  ovlap[, c("chr", "ori", "pos1", "seq1") := rtea[idx1, .(chr, ori, pos, seq)]]
  ovlap[, c("pos2", "seq2") := rtea[idx2, .(pos, seq)]]
  ovlap[, overhangint1 := overhangint[idx1]]
  ovlap[, overhangint2 := overhangint[idx2]]
  ovlap[, ungappos1 := ungappos[idx1]]
  ovlap[, ungappos2 := ungappos[idx2]]
  ovlap %<>% merge(baseseqdt,
                   by.x = c("chr", "pos1", "ori", "overhangint1", "ungappos1"),
                   by.y = c("chr", "pos", "ori", "overhangint", "ungappos"))
  setnames(ovlap, "seq", "baseseq1")
  ovlap %<>% merge(baseseqdt,
                   by.x = c("chr", "pos2", "ori", "overhangint2", "ungappos2"),
                   by.y = c("chr", "pos", "ori", "overhangint", "ungappos"))
  setnames(ovlap, "seq", "baseseq2")
  ovlap[, c("chr", "pos1", "pos2", "overhangint1", "overhangint2", "ungappos1", "ungappos2") := NULL]
  ovlap[ori == "f" & posgap > 0, seq1 := paste0(seq1, stringi::stri_sub(baseseq1, 1, posgap))]
  ovlap[ori == "f" & posgap < 0, seq2 := paste0(seq2, stringi::stri_sub(baseseq2, 1, posgap))]
  ovlap[ori == "r" & posgap > 0, seq2 := paste0(stringi::stri_sub(baseseq2, -posgap), seq2)]
  ovlap[ori == "r" & posgap < 0, seq1 := paste0(stringi::stri_sub(baseseq1, -posgap), seq1)]
  ovlap[, c("baseseq1", "baseseq2") := NULL]
  ovlap %<>% .[, .(idx1 = list(idx1), idx2 = list(idx2)), by = .(seq1, seq2, ori)]
  # comparing sequences
  if(verbose) msg("Comparing sequences...")
  ovlap$seqmatch <- bptry(bpvec(
    seq_len(nrow(ovlap)),
    function(i) {
      gc()
      ovlap[i, pSeqMatchRatio(seq1, seq2, ori, shift_range = seqmatch_shift_range)]
    },
    BPPARAM = MulticoreParam(workers = threads, tasks = nrow(ovlap) / 1e5, log = verbose)
  ))
  ovlap %<>% .[seqmatch >= seqmatch_cutoff]
  
  # finding lowest idx
  if(verbose) msg("Finding lowest idx...")
  ovlist <- ovlap[, S4Vectors::SelfHits(c(seq_len(nrow(rtea)), unlist(idx1), unlist(idx2)),
                                        c(seq_len(nrow(rtea)), unlist(idx2), unlist(idx1)),
                                        nrow(rtea))] %>% as.list
  lowgrp <- mclapply(ovlist, function(x) min(unlist(ovlist[x]))) %>% unlist
  grp <- mclapply(ovlist, function(x) min(lowgrp[x])) %>% unlist
  while(!identical(lowgrp, grp)) {
    lowgrp <- grp
    grp <- mclapply(ovlist, function(x) min(lowgrp[x])) %>% unlist
  }
  dupgrp <- which(table(grp) > 1) %>% names %>% as.integer
  
  if(selectbest) {
    if(verbose) msg("Selecting best idx...")
    bestone <- mclapply(dupgrp, function(g) {
      idx <- which(grp == g)
      passed <- idx[rtea[idx, Filter == "PASS"]]
      if(length(passed) == 0) passed <- idx
      highscore <- rtea[passed, TEscore] %>% {
        passed[. == suppressWarnings(max(., na.rm = T))]} %>%
        na.omit
      if(length(highscore) == 0) highscore <- passed
      highcnt <- rtea[highscore, uniqueCnt] %>% {highscore[. == max(.)]}
      prox <- if(rtea[highcnt, ori][1] == "f") {
        rtea[highcnt, pos] %>% {highcnt[. == max(.)]}
      } else {
        rtea[highcnt, pos] %>% {highcnt[. == min(.)]}
      }
      prox[1]
    }) %>% unlist
    bestone[match(grp, dupgrp)]
  } else {
    grp[!grp %in% dupgrp] <- NA
    grp
  }
  
}

matchScallop.ctea <- function(ctea, scallopfile, matchrange = 300, overhangmin = 5) {
  require(rtracklayer)
  scallop <- import(scallopfile, "gtf")
  strand(scallop) <- "*"
  trspt <- subset(scallop, type == "transcript")
  exon <- subset(scallop, type == "exon")
  m <- match(exon$transcript_id, trspt$transcript_id)
  exon$trstart <- start(trspt[m])
  exon$trend <- end(trspt[m])
  firstexon <- subset(exon, start == trstart, select = transcript_id)
  fgr <- shift(firstexon, -1) %>%
    resize(width = pmin(width(.), matchrange))
  lastexon <- subset(exon, end == trend, select = transcript_id)
  rgr <- shift(lastexon, 1) %>%
    resize(width = pmin(width(.), matchrange), fix = "end")
  ctea[ori == "f", scallop_id := grLocator(chr, pos, refgr = fgr)]
  ctea[ori == "f" & overhang <= overhangmin,
       scallop_id := grLocator(chr, pos + gap, refgr = fgr)]
  ctea[ori == "r", scallop_id := grLocator(chr, pos, refgr = rgr)]
  ctea[ori == "r" & overhang <= overhangmin,
       scallop_id := grLocator(chr, pos - gap, refgr = rgr)]
  ctea
}

annotateScallop.ctea <- function(rtea, scallopfile, edbpkg = edbPKG, threads = getOption("mc.cores", 2)) {
  require(rtracklayer)

  msg("Loading exons...")
  library(edbpkg, character.only = T)
  edb <- get(edbpkg)
  exons <- exons(edb, columns = "tx_id")
  exons %<>% .[order(.$tx_id)]
  msg("Annotating scallop transcript...")

  scallop <- import(scallopfile, "gtf")
  seqlevels(scallop) %<>% sub("chr", "", .)
  noNA <- rtea[, which(!is.na(scallop_id))]
  n <- nrow(rtea)
  tx_id <- mclapply(
    noNA,
    mc.cores = threads,
    function(idx) {
      cat(sprintf("\r%0.2f%%              ", idx/n*100))

      scallop_id <- strsplit(rtea[idx, scallop_id], ",")[[1]]
      clippos <- ungapPos.rtea(rtea[idx])
      ori <- rtea[idx, ori]
      hardpos <- rtea[idx, ifelse(ori == "f", hardstart, hardend)]
      tepos <- ifelse(is.na(hardpos), clippos, hardpos)

      whichTranscript <- function(scl, overlapprop_cutoff = 0.7) {
        scl %<>% subset(type == "exon")
        ovtid <- subsetByOverlaps(exons, scl)$tx_id %>% unique
        if(length(ovtid) == 0) return(NA)
        overlapsize <- sapply(ovtid,
                              function(x) subset(exons, tx_id == x) %>%
                                intersect(scl) %>%
                                width %>%
                                sum
        )
        bestid <- overlapsize %>%
          .[. > sum(width(scl)) * overlapprop_cutoff] %>%
          sort(decreasing = T) %>%
          names %>%
          .[1]
      }
      whichTranscript(subset(scallop, transcript_id %in% scallop_id))
  }) %>% unlist
  ft <- fusiontype(rtea[noNA], tx_id)[, !names(rtea), with = F]
  names(ft) %<>% paste0("scallop_", .)
  idx <- rep(NA, n)
  idx[noNA] <- seq_along(noNA)
  data.table(rtea, ft[idx])
}

# hardstart, hardend column is needed on rtea. If not present run localHardClip.
fusiontype <- function(rtea, tx_id = rtea$fusion_tx_id,
                       edbpkg = edbPKG,
                       columns = c("tx_biotype", "gene_id", "gene_name")) {
  stopifnot(nrow(rtea) == length(tx_id))

  require(edbpkg, character.only = T)
  edb <- get(edbpkg)

  tx <- if(sum(!is.na(tx_id)) > 0) {
    transcripts(edb,
                columns = c("seq_strand", "tx_seq_start", "tx_seq_end", columns),
                filter = TxIdFilter(na.omit(tx_id)),
                return.type = "GRanges") %>% data.frame %>% data.table
    # return.type = "data.frame") %>% data.table
  } else {
    data.table(tx_id = NA, start = NA, end = NA, strand = NA,
               as.list(rep(NA, length(columns))) %>% structure(names = columns) %>% do.call(data.table, .))
  }
  dt <- data.table(ugpos = ungapPos.rtea(rtea),
                   rtea[, .(ori, hardstart, hardend)],
                   tx[tx_id, on = "tx_id"]
  )
  dt[is.na(tx_id), fusion_type := "intergenic"]
  dt[(ori == "f" & ugpos > end) | (ori == "r" & ugpos < start), fusion_type := "impossible"]
  dt[between(ugpos, start, end, NAbounds = NA), fusion_type := "exonic/exonization"]
  dt[(ori == "f" & pmin(ugpos, hardstart, na.rm = T) < start),
     fusion_type := fifelse(strand == "+", "alternative TSS", "read-through")]
  dt[(ori == "r" & pmax(ugpos, hardend, na.rm = T) > end),
     fusion_type := fifelse(strand == "+", "read-through", "alternative TSS")]
  dt <- dt[, c("fusion_type", "tx_id", columns), with = F]
  names(dt)[-1] %<>% paste0("fusion_", .)
  if(exists("fusion_tx_id", rtea) && identical(rtea$fusion_tx_id, dt$fusion_tx_id)) {
    dt$fusion_tx_id <- NULL
  }
  data.table(rtea, dt)
}

fusionTxMatch <- function(sam, edbpkg, strandedness = c("non-stranded", "first-sense", "first-antisense")) {
  if(length(sam) == 0) {
    return(data.table(tx_id = NA_character_,
                      tx_support_exon = NA_integer_,
                      tx_support_intron = NA_integer_,
                      numgap = NA_integer_)
    )
  }

  library(edbpkg, character.only = T)
  edb <- get(edbpkg)

  splvec <- function(x, FUN, maxN, AGGREGATE = c, ...) {
    spl <- lapply(seq(1, length(x), maxN), function(i1) i1:min(i1 + maxN - 1, length(x)))
    lst <- lapply(spl, function(i) FUN(x[i], ...))
    do.call(AGGREGATE, lst)
  }
  maxlen <- 900

  if(missing(strandedness)) {
    strandedness <- "non-stranded"
  }
  rnastrand <- function(gr, flag, stranded = strandedness) {
    if(stranded == "non-stranded") {
      unstrand(gr)
    } else if(stranded == "first-antisense") {
      strand(gr)[bamFlagTest(flag, "isFirstMateRead")] %<>% invertStrand
      gr
    } else if(stranded == "first-sense") {
      strand(gr)[bamFlagTest(flag, "isSecondMateRead")] %<>% invertStrand
      gr
    } else {
      stop("The value of stranded should be 'non-stranded', 'first-sense', or 'first-antisense', but is: ", stranded)
    }
  }

  irl <- extractAlignmentRangesOnReference(cigar(sam), pos = start(sam))
  gr <- GRanges(rep(seqnames(sam), elementNROWS(irl)),
                unlist(irl),
                rep(strand(sam), elementNROWS(irl))
  )
  gr %<>% rnastrand(rep(mcols(sam)$flag, elementNROWS(irl)))
  exons <- if(length(gr) <= maxlen) {
    exons(edb, columns = "tx_id", filter = GRangesFilter(reduce(gr)))
  } else {
    exons <- splvec(reduce(gr),
                    function(x) exons(edb, columns = "tx_id", filter = GRangesFilter(x)),
                    maxN = maxlen
    )
    exons[!duplicated(data.frame(exons))]
  }

  ovlexon <- data.table(
    tx_id = exons$tx_id,
    cnt = findOverlaps(narrow(gr[width(gr) > 7], 5, -5), exons, type = "within") %>% countSubjectHits
  ) %>% .[, .(tx_support_exon = sum(cnt)), by = tx_id]

  gap <- junctions(sam) %>% rnastrand(mcols(sam)$flag) %>% unlist
  gaplen <- length(unique(gap))
  introns <- if(gaplen > 0) {
    if(gaplen <= maxlen) {
      intronsByTranscript(edb, filter = GRangesFilter(reduce(gap))) %>% unlist
    } else {
      introns <- splvec(reduce(gap),
                        function(x) intronsByTranscript(edb, filter = GRangesFilter(x)),
                        maxN = maxlen
      ) %>% unlist
      introns[!duplicated(data.frame(introns))]
    }
  } else {
    GRanges()
  }
  ovlintron <- if(length(introns) > 0) {
    data.table(
      tx_id = names(introns),
      cnt = countOverlaps(introns, gap, type = "equal")
    ) %>% .[, .(tx_support_intron = sum(cnt)), by = tx_id]
  } else {
    data.table(tx_id = character(0), tx_support_intron = integer(0))
  }
  ovlcnt <- merge(ovlexon, ovlintron, all = T, by = "tx_id")
  ovlcnt[is.na(tx_support_exon), tx_support_exon := 0]
  ovlcnt[is.na(tx_support_intron), tx_support_intron := 0]
  ovlcnt %<>% .[tx_support_exon > 0 | tx_support_intron > 0]
  data.table(ovlcnt[order(-tx_support_intron, -tx_support_exon, tx_id)][1],
             numgap = length(gap))
}

fusiontypeByCigar <- function(rtea, bamfile,
                              strandedness = c("non-stranded", "first-sense", "first-antisense"),
                              edbpkg = edbPKG,
                              threads = getOption("mc.cores", detectCores())
) {
  require(BiocParallel)
  require(GenomicAlignments)
  require(ensembldb)

  if(missing(strandedness)) {
    strandedness <- checkStrandedness(bamfile, edbpkg)
  }

  if(paste0("package:", edbpkg) %in% search()) {
    detach(paste0("package:", edbpkg), unload = T, character.only = T)
  }

  trptMatch <- function(i, edbpkg, bamfile, chr, pos, ori) {
    cat(sprintf("\r%0.2f%%              ", i/n*100))
    sam <- getClippedPairs(bamfile, chr, pos, ori, subsample = T)
    sam <- c(first(sam$pairs), second(sam$pairs), sam$nopair)
    fusionTxMatch(sam, edbpkg, strandedness)
  }
  msg("Finding matching transcripts...")
  n <- nrow(rtea)
  trpt <- bptry(
    bpmapply(
      FUN = trptMatch,
      seq_len(n),
      rtea$chr, rtea$pos, rtea$ori,
      MoreArgs = list(edbpkg = edbpkg, bamfile = bamfile),
      SIMPLIFY = F,
      BPPARAM = MulticoreParam(workers = threads, stop.on.error = T, tasks = n)
    ),
    bplist_error = function(e) {
      message(e, "\n")
      save(rtea, bamfile, threads, e,
           file = "fusiontypeByJunctionErr.RData")
      e
    }
  ) %>% rbindlist

  rtea <- fusiontype(rtea, trpt$tx_id, edbpkg = edbpkg)
  stopifnot(identical(trpt$tx_id, rtea$fusion_tx_id))
  data.table(rtea, trpt[, .(tx_support_exon, tx_support_intron, numgap)])
}

nearbypair <- function(rtea,
                          overhang_cutoff = 10L,
                          maxTSD = 10,
                          maxgap = 2e5L) {
  library(GenomicRanges)

  rtea[, ungappos := ungapPos.rtea(rtea, overhang_cutoff = overhang_cutoff)]
  rgrr <- rtea[ori == "r", GRanges(chr, IRanges(ungappos, ungappos))]
  rgrf <- rtea[ori == "f", GRanges(chr, IRanges(ungappos, ungappos))]
  ovlap <- findOverlaps(rgrr, rgrf, maxgap = maxgap, select = "all")
  ovlap %<>% .[rtea[ori == "r"][queryHits(.), ungappos] < rtea[ori == "f"][subjectHits(.), ungappos] + maxTSD]
  pairr <- rtea[ori == "r"][queryHits(ovlap)]
  pairf <- rtea[ori == "f"][subjectHits(ovlap)]
  meta <- data.table(
    dist = pairf$ungappos - pairr$ungappos,
    sameGene = pairr$gene_name == pairf$gene_name,
    sameSource = abs(pairf$hardstart - pairr$hardend) < 6000 & pairf$hardTE == pairf$hardTE,
    classr = pairr$class,
    classf = pairf$class
  )
  meta[pairr$class == pairf$class, class := classr]
  meta[pairr$class != pairf$class, class := "different"]
  meta[pairr$isPolyA == T, class := classf]
  meta[pairf$isPolyA == T & classr != "PolyA", class := classr]
  meta[pairr$isPolyA & !pairf$isPolyA, polyA := "r"]
  meta[!pairr$isPolyA & pairf$isPolyA, polyA := "f"]
  meta[pairr$isPolyA & pairf$isPolyA, polyA := "both"]
  meta[!pairr$isPolyA & !pairf$isPolyA, polyA := "none"]
  meta[, Passed := F]
  meta[pairr$Filter == "PASS" & pairf$Filter == "PASS", Passed := T]
  meta[pairr$Filter == "PASS" & pairf$isPolyA, Passed := T]
  meta[pairf$Filter == "PASS" & pairr$isPolyA, Passed := T]
  meta[, classr := NULL]
  meta[, classf := NULL]

  list(r = pairr,
       f = pairf,
       meta = meta
  )

}

geneticLocation <- function(chr, start, end = start) {
  # require(rtracklayer)
  # writeLines(paste("Importing", gene_file))
  # annodata <- import(gene_file, "gtf")
  annodata <- readRDS(gene_file)
  msg("Location idenfying")
  txannodata <- subset(annodata, type != "gene")
  stopifnot(length(chr) == length(start))
  isNA <- is.na(chr) | is.na(start) | is.na(end)
  ch <- chr[!isNA]
  st <- pmin(start, end)[!isNA]
  en <- pmax(start, end)[!isNA]
  chr <- sub("chr", "", chr)
  cgr <- GRanges(ch, IRanges(st, en), "*")
  overlaptx <- findOverlaps(cgr, txannodata, ignore.strand = T)
  selectAnnodata <- function(overlap, data) {
    matchdata <- data.table(idx = queryHits(overlap),
                            data[subjectHits(overlap)] %>%
                              data.frame
    )
    matchdata[, tx_order := sub(".*-", "", transcript_name) %>% as.integer]
    matchdata <- matchdata[
      ,
      .SD[tx_order == min(tx_order)],
      by = idx,
      .SDcols = strand:exon_number
      ]
    matchdata[, gene_order := substring(gene_id, 5) %>% as.integer]
    matchdata[
      ,
      .SD[which.min(gene_order)],
      by = idx,
      .SDcols = strand:exon_number
      ]
  }
  txmatch <- selectAnnodata(overlaptx, txannodata)
  anno <- txmatch[match(seq_len(sum(!isNA)), idx), !"idx"]
  anno[is.na(type), type := "intergenic"]
  anno[!is.na(exon_number), type_number := paste(type, exon_number, sep="_")]
  allidx <- rep(NA_integer_, length(start))
  allidx[!isNA] <- seq_len(sum(!isNA))
  anno[allidx, .(gene_id, gene_name, transcript_id, transcript_name, type, type_number, strand)]
}

annotate.ctea <- function(ctea,
                          flank_width = 100000L,
                          junction_width = 5L,
                          overhang_cutoff = junction_width
                 ) {
  # require(rtracklayer)
  # writeLines(paste("Importing", gene_file))
  # annodata <- import(gene_file, "gtf")
  library(GenomicRanges)
  annodata <- readRDS(gene_file)
  seqlevels(annodata) %<>% pastechr
  annodata$exon_number <- as.integer(annodata$exon_number)
  msg("Annotating TEI.")
  txannodata <- subset(annodata, type != "gene")
  ungappos <- ungapPos.rtea(ctea, overhang_cutoff)
  cgr <- ctea[, GRanges(pastechr(chr), IRanges(ungappos, width=1), c(r = "+", f = "-")[ori])]
  overlaptx <- findOverlaps(cgr, txannodata, ignore.strand = T)
  selectAnnodata <- function(overlap, data) {
    matchdata <- data.table(idx = queryHits(overlap),
                            data[subjectHits(overlap)] %>% data.frame
                 )
    matchdata[, tx_order := sub(".*-", "", transcript_name) %>% as.integer]
    matchdata <- matchdata[,
                           .SD[tx_order == min(tx_order)],
                           by = idx,
                           .SDcols = strand:exon_number
                 ]
    matchdata[, gene_order := substring(gene_id, 5) %>% as.integer]
    matchdata[,
              .SD[which.min(gene_order)],
              by = idx,
              .SDcols = strand:exon_number
    ]
  }
  txmatch <- selectAnnodata(overlaptx, txannodata)
  anno <- txmatch[match(seq_len(nrow(ctea)), idx), !"idx"]
  anno[is.na(type), type := "intergenic"]
  rm(txannodata)

  spldonor <- subset(annodata, type == "intron") %>%
    flank(width = junction_width, start = T, both = T)
  mcols(spldonor)$type <- "splice_donor"
  overlapdo <- findOverlaps(cgr, spldonor, ignore.strand = F)
  domatch <- selectAnnodata(overlapdo, spldonor)
  anno[domatch$idx] <- domatch[, !"idx"]
  rm(spldonor)

  splacceptor <- subset(annodata, type == "intron") %>%
    flank(width = junction_width, start = F, both = T)
  mcols(splacceptor)$type <- "splice_acceptor"
  levels(strand(cgr)) <- c("+" = "-", "-" = "+", "*" = "*")[levels(strand(cgr))]
  overlapac <- findOverlaps(cgr, splacceptor, ignore.strand = F)
  acmatch <- selectAnnodata(overlapac, splacceptor)
  anno[acmatch$idx] <- acmatch[, !"idx"]
  rm(splacceptor)

  anno[!is.na(exon_number), type_number := paste(type, exon_number, sep="_")]
  anno[, geneSide := ifelse(xor(ctea$ori == "f", strand == "+"), "5", "3")]

  upstream <- flank(subset(annodata, type == "gene"),
                    width = flank_width,
                    start = T
  )
  mcols(upstream)$type <- "upstream"
  overlapup <- findOverlaps(cgr, upstream, ignore.strand = T)
  upmatch <- data.table(idx = queryHits(overlapup),
                        gene = upstream[subjectHits(overlapup)] %>%
                          mcols %>%
                          .$gene_name
  )
  upmatch <- upmatch[, .(gene = paste(gene, collapse = ",")), by = idx]
  anno[, upstream := upmatch[match(.I, idx), gene]]

  downstream <- flank(subset(annodata, type == "gene"),
                      width = flank_width,
                      start = F
  )
  mcols(downstream)$type <- "downstream"
  overlapdw <- findOverlaps(cgr, downstream, ignore.strand = T)
  dwmatch <- data.table(idx = queryHits(overlapdw),
                        gene = upstream[subjectHits(overlapdw)] %>%
                          mcols %>%
                          .$gene_name
  )
  dwmatch <- dwmatch[, .(gene = paste(gene, collapse = ",")), by = idx]
  anno[, downstream := dwmatch[match(.I, idx), gene]]
  names(anno) %<>% paste0("pos_", .)
  data.table(ctea, anno[, .(pos_gene_id, pos_gene_name, pos_transcript_id, pos_transcript_name,
                            pos_type, pos_type_number, pos_strand, pos_upstream, pos_downstream)
                        ]
  )

}

localHardClip <- function(rtea,
                          mateDistMax = 100000,
                          overhangmin = 5,
                          score_cutoff = 10,
                          threads = getOption("mc.cores", detectCores())) {
  if(nrow(rtea) == 0) {
    salign <- data.table(hardstart = integer(0),
                         hardend = integer(0),
                         hardTE = character(0),
                         hardDist = integer(0),
                         hardNumIntron = integer(0),
                         hardSpl = character(0)
    )
    return( data.table(rtea, salign) )
  }
  require(bsgenomePKG, character.only = T)
  genome <- get(bsgenomePKG)
  searchstart <- searchend <- rtea$pos  # ungapPos.rtea(rtea)
  searchend[which(rtea$overhang < overhangmin & rtea$ori == "f")] <- rtea[overhang < overhangmin & ori == "f", pos + gap]
  searchstart[which(rtea$overhang < overhangmin & rtea$ori == "r")] <- rtea[overhang < overhangmin & ori == "r", pos - gap]
  searchstart[rtea$ori == "f"] <- rtea[ori == "f", pmax(pos - mateDist - 100, pos - mateDistMax, 1)]
  searchend[rtea$ori == "r"] <- rtea[ori == "r", pos + pmin(mateDist + 100, mateDistMax)]
  searchend <- pmin(searchend, seqlengths(genome)[pastechr(rtea$chr)])
  refseq <- getSeq(genome, pastechr(rtea$chr), searchstart, searchend)
  refseq %<>% DNAStringSet
  matchscore <- nucleotideSubstitutionMatrix(match = 1, mismatch = -10, baseOnly = TRUE) %>%
    cbind(N = 0) %>%
    rbind(N = 0)
  n <- rtea[, .N]
  msg("Hard clip align ...")
  salign <- mclapply(
    seq_len(n),
    mc.cores = threads,
    mc.preschedule = T,
    function(i) {
      cat(sprintf("\r%0.2f%%              ", i/n*100))
      rs <- refseq[i]
      cs <- rtea[i, seq]
      align <- pairwiseAlignment(cs, rs, substitutionMatrix = matchscore, type = "global-local")
      list(score = score(align),
           mstart = start(subject(align)),
           mend = end(subject(align)))
    }
  ) %>% rbindlist

  msg("Done align.")
  salign[, chr := rtea$chr]
  salign[, hardstart := searchstart + mstart - 1]
  salign[, hardend := searchstart + mend - 1]
  salign[score < score_cutoff, hardstart := NA]
  salign[score < score_cutoff, hardend := NA]
  reploc <- salign[score >= score_cutoff, repeatLocation(chr, hardstart, hardend)]
  salign[, hardTE := NA_character_]
  salign[score >= score_cutoff, hardTE := reploc[[2]]]
  salign[, hardDist := ifelse(rtea$ori == "f", searchend - hardend, hardstart - searchstart)]

  # require(rtracklayer)
  # writeLines(paste("Importing", gene_file))
  # annodata <- import(gene_file, "gtf")

  hclipped <- salign$score >= score_cutoff
  # number of exons between mapping position and hard clip
  if(any(hclipped)){
    annodata <- readRDS(gene_file)
    seqlevels(annodata) %<>% pastechr
    st <- searchstart
    en <- searchend
    st[rtea$ori == "f"] <- salign$hardstart[rtea$ori == "f"]
    en[rtea$ori == "r"] <- salign$hardend[rtea$ori == "r"]
    gr <- GRanges(pastechr(rtea$chr[hclipped]), IRanges(st[hclipped], en[hclipped]), "*")
    overlap <- findOverlaps(gr, annodata, ignore.strand = T)
    if(length(overlap) == 0) {
      salign[, hardNumIntron := NA_integer_]
    } else {
      overlapdt <- data.table(idx = queryHits(overlap),
                              type = annodata$type[subjectHits(overlap)],
                              transcript = annodata$transcript_id[subjectHits(overlap)])
      whichenst <- overlapdt[, .N, .(idx, transcript)][order(idx, -N, transcript)][, .SD[1], by = idx]
      whichenst <- whichenst[match(seq_len(max(whichenst$idx)), idx), transcript]
      numint <- overlapdt[, .(num_intron = sum(type[transcript == whichenst[idx]] == "exon", na.rm = T)), by = idx]
      salign$hardNumIntron <- numint[match(seq_len(n), which(hclipped)[idx]), num_intron]
      salign[is.na(hardNumIntron) & !is.na(hardstart), hardNumIntron := 0]
    }

  } else {
    salign[, hardNumIntron := NA_integer_]
  }

  # whether the hard clip is on the normal splicingn site
  splclipped <- hclipped & rtea$type != "intergenic"
  if(any(splclipped)) {
    ungappos <- ungapPos.rtea(rtea[splclipped == T])
    gr <- GRanges(pastechr(rtea[splclipped == T, chr]), IRanges(ungappos, ungappos), "*")
    intrdata <- subset(annodata, type == "intron")
    spldata <- c(flank(intrdata, width = 5L, start = T, both = T),
                 flank(intrdata, width = 5L, start = F, both = T)
                 )
    intovlap <- findOverlaps(gr, spldata, maxgap = 5L) %>%
      data.frame %>%
      data.table
    intovlap[, hardDist := salign[splclipped == T][queryHits, hardDist]]
    intovlap[, width := rep(width(intrdata), 2)[subjectHits]]
    intovlap[, enst := spldata$transcript_id[subjectHits]]
    intovlap[, rdiff := abs(hardDist - width)]
    setorder(intovlap, enst)
    hardmatched <- intovlap[rdiff < 10L,
                            .(enst = enst[which.min(rdiff)]),
                            by = queryHits
                            ]
    salign$hardSpl <- hardmatched[match(seq_len(n), which(splclipped)[queryHits]), enst]
  } else {
    salign[, hardSpl := NA_character_]
  }

  data.table(rtea, salign[, hardstart:hardSpl])
}

TEcoordinate <- function(rtea, threads = getOption("mc.cores", detectCores())) {
  require(Biostrings)
  require(parallel)
  msg("Reading TE fasta file...")
  fa <- readDNAStringSet(refTEfa)
  names(fa) %<>% sub("\\s.*", "", .)
  TEclass <- rtea[, family]
  TEclass[TEclass == "L1HS"] <- "LINE1"
  mfa <- match(TEclass, names(fa))
  idxfa <- as.list(mfa)
  idxfa[TEclass == "LINE1"] <- list(grep("L1", names(fa)))
  idxfa[TEclass == "AluY"] <- list(grep("ALU", toupper(names(fa))))
  idxfa[TEclass == "SVA"] <- list(grep("SVA", names(fa)))
  # idxfa[TEclass == "HERV"] <- list(grep("HERV|MER|LTR|MLT", names(fa)))
  idxfa[TEclass %in% c("HERVK", "LTR5_Hs", "LTR5A", "LTR5B", "LTR5")] %<>% lapply(c, match("HERVK-full", names(fa)))
  idxfa[TEclass %in% c("HERVH", "LTR7", "LTR7Y", "LTR7B", "LTR7C")] %<>% lapply(c, match("HERVH-full", names(fa)))
  isna <- is.na(sapply(idxfa, `[`, 1))
  idxfa <- idxfa[!isna]
  fseq <- rtea[!isna, DNAStringSet(seq)]
  rseq <- reverseComplement(fseq)
  ori <- rtea[!isna, ori]
  msg("Mapping + strand to TE references...")
  n <- sum(!isna)
  falign <- mclapply(
    seq_len(n),
    mc.cores = threads,
    mc.preschedule = T,
    function(i) {
      cat(sprintf("\r%0.2f%%              ", i/n*100))
      pairwiseAlignment(fa[idxfa[[i]]],
                        fseq[[i]],
                        type = "local-global",
                        gapExtension = 2
      )
    }
  )

  msg("Mapping - strand to TE references ...")
  ralign <- mclapply(
    seq_len(n),
    mc.cores = threads,
    mc.preschedule = T,
    function(i) {
      cat(sprintf("\r%0.2f%%              ", i/n*100))
      pairwiseAlignment(fa[idxfa[[i]]],
                        rseq[[i]],
                        type = "local-global",
                        gapExtension = 2
      )
    }
  )
  msg("Mapping done.")
  fscore <- lapply(falign, score)
  rscore <- lapply(ralign, score)
  fscomax <- sapply(fscore, max)
  rscomax <- sapply(rscore, max)
  isforward <- fscomax >= rscomax
  TEscore <- ifelse(isforward, fscomax, rscomax)
  TEside <- rep("5", n)
  TEside[xor(isforward, ori == "f")] <- "3"
  TEbreak <- list()
  TEbreak[TEside == "5"] <- lapply(ifelse(isforward, falign, ralign)[TEside == "5"],
                                   function(x) end(pattern(x))
  )
  TEbreak[TEside == "3"] <- lapply(ifelse(isforward, falign, ralign)[TEside == "3"],
                                   function(x) start(pattern(x))
  )
  bestmatch <- list()
  bestmatch[isforward] <- Map(function(x, y) which(x==y),
                              fscore[isforward],
                              fscomax[isforward]
                          )
  bestmatch[!isforward] <- Map(function(x, y) which(x==y),
                               rscore[!isforward],
                               rscomax[!isforward]
                           )
  TEbreak <- mapply(function(x, i) paste(x[i], collapse = ","), TEbreak, bestmatch)
  TEfamily <- mapply(function(i, j) paste(names(fa[i[j]]), collapse = ","), idxfa, bestmatch)
  idx <- rep(NA_integer_, nrow(rtea))
  idx[!isna] <- seq_len(n)
  data.table(rtea, data.table(TEfamily, TEscore, TEside, TEbreak)[idx, ])
}

writeIGVbat <- function(bams,
                        chr, pos, width = 150,
                        outdir = getwd(),
                        bat = "igv.bat", append = T) {
  rmskfile <- "/home/bl177/ref/TEtranscripts/GRCh37_rmsk.gtf"
  stopifnot(length(chr)==length(pos))
  stopifnot(all(file.exists(bams)))
  start <- pos - round(width / 2)
  end <- pos + round(width / 2)
  outFile <- file(bat, ifelse(append, "a", "w"))
  writeLines("new", outFile)
  writeLines(paste("load", rmskfile), outFile)
  writeLines(paste("load", bams), outFile)
  writeLines(paste("snapshotDirectory", outdir), outFile)
  for(i in seq_along(chr)){
    writeLines(paste0("goto ", chr[i], ":", start[i], "-", end[i]), outFile)
    writeLines("sort base", outFile)
    prefix <- sub(".bam$", "", basename(bams[1]))
    pngFile <- paste0(prefix, ".",
                      chr[i], "_",
                      pos[i], ".png")
    writeLines(paste("snapshot", pngFile), outFile)
  }
  close(outFile)
  invisible(0)
}

IGVsnapshot.rnatea <- function(rnatea, outdir, bamfiles, wgsbamfiles = NULL) {
  bamfiles[] %<>% path.expand
  if(!missing(wgsbamfiles)) wgsbamfiles[] %<>% path.expand
  if(!file.exists(outdir)) dir.create(outdir, recursive = T)
  igvbat <- paste0(outdir, "/igv_", Sys.Date(), ".bat")
  if(file.exists(igvbat)) {
    unlink(igvbat)
  }
  rnatea[, writeIGVbat(na.omit(c(bamfiles[ID], wgsbamfiles[ID])),
                       chr, pos,
                       outdir = outdir, bat = igvbat), by = ID]
  con <- file(igvbat, "a")
  writeLines("exit", con)
  close(con)
  # system(paste("java -jar $IGV/igv.jar -b", igvbat))
  system(paste("xvfb-run igv -b", igvbat))
}

readFastp <- function(fastpfiles) {
  require(data.table)
  require(rjson)
  lapply(fastpfiles, function(x) {
    qclist <- fromJSON(file = x)
    qc <- qclist$summary$before_filtering
    qc %<>% {c(.,
               filtered_reads = .[["total_reads"]] - qclist$summary$after_filtering$total_reads,
               dup_rate = qclist$duplication$rate,
               trimmed_reads = qclist$adapter_cutting$adapter_trimmed_reads,
               trimmed_bases = qclist$adapter_cutting$adapter_trimmed_bases,
               insert_size = qclist$insert_size$peak
    )}
  }) %>% rbindlist
}
