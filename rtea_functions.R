# 190624 - not filtering inserted polyA tail (filterSimpleRepeat.ctea function)
# 190814 - duplicate removal by read name
# 190829 - TEcoordinate reference fasta update (L1, L1HS)
# 190909 - considering R1 or R2 in duplicate removal
# 191121 - considering only clipped length and mate position (not clipped sequence) in marking possible duplicate
# 191121 - repeat location annotation bug fix (when family is NA then class is assinged)
# 191127 - filter out normal polyA site
# 191127 - flag low TEscore (cntFilter.ctea function)
# 191129 - considering same clipped seq also as possible duplicate
# 191204 - not filtering polyA when reading ctea file (readctea function)
# 191223 - converted gtf files to rds files to reduce file size.
# 200103 - before removing duplicate reads, filter out the reads of which start(ori == f) or end(ori == r) site is too far away.

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
# polyTE_file <- paste0("/home/bl177/lee/boram/ref/polymorphicTE/GRanges_union_non_reference_TEI_", build, ".rds")
polyTE_file <- file.path(refdir, "GRanges_union_non_reference_TEI.rds")
bsgenomePKG <- paste0("BSgenome.Hsapiens.UCSC.", build)
# gene_file <- file.path(refdir, "annotation_data.gtf")
gene_file <- file.path(refdir, "annotation_data.rds")
# rmsk_gtf_file <- file.path(refdir, "rmsk.gtf")
rmsk_file <- file.path(refdir, "rmsk.rds")
refTEfa <- getOption("refTEfa", file.path(refdir, "consensus_L1_ALU_SVA_HERV.fa"))
# if(is.null(refTEfa)) refTEfa <- "/home/bl177/rnatea/ref/hg38.RepeatMasker-4.0.6-Dfam-2.0.fa"

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
  ctea[substr(family, 1, 4) == "HERV", class := "HERV"]

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
                             Alength = 10, percentA = 0.8,
                             refAlength = 20, refPercentA = 0.8) {
  if(!exists("posRep", ctea)) {
    ctea <- repeatPosition.ctea(ctea)
  }
  require(bsgenomePKG,  character.only = T)
  genome <- get(bsgenomePKG)
  st <- ctea[, pmax(ifelse(ori == "f", pos + 1, pos - refAlength), 
                    0
                    )
             ]
  en <- ctea[, pmin(ifelse(ori == "f", pos + refAlength, pos - 1), 
                    seqlengths(genome)[pastechr(chr)]
                    )
             ]
  refrootseq <- ctea[, getSeq(genome, pastechr(chr), st, en)]
  
  Asite <- stringr::str_count(refrootseq, "A") / refAlength >= refPercentA
  Tsite <- stringr::str_count(refrootseq, "T") / refAlength >= refPercentA
  
  proxseq <- ctea[, ifelse(ori == "f", 
                           stringr::str_sub(seq, start = -Alength),
                           stringr::str_sub(seq, end = Alength)
  )]
  Aprop <- stringr::str_count(proxseq, "A") / Alength
  Tprop <- stringr::str_count(proxseq, "T") / Alength
  simpleSite <- (Asite & Aprop >= percentA) | (Tsite & Tprop >= percentA)
  ctea[simpleSite == F]
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

clippedBam <- function(bamfile, gr, mapqFilter = 1L, yieldSize = 1e5) {
  library(GenomicAlignments)
  library(GenomicFiles)
  
  cbfile <- filterBam(bamfile, tempfile(),
                      param = ScanBamParam(which = gr, mapqFilter = mapqFilter),
                      indexDestination = F
  )
  on.exit(unlink(cbfile))
  
  bf <- BamFile(cbfile, yieldSize = yieldSize)
  yield <- function(x)
    readGAlignments(x,
                    param = ScanBamParam(
                      what = c("qname", "seq", "qual", "mapq", "flag", "mpos"),
                      tag = "NM",
                      mapqFilter = mapqFilter
                    )
    )
  reduceByYield(bf, yield, identity, REDUCEsampler(yieldSize, F))
}

getClippedReads <- function(bamfile, chr, pos, ori = c("f", "r"), 
                            searchWidth = 10L,
                            mapqFilter = 1L,
                            subsample = F,
                            maxReads = 1e5) {
  require(GenomicAlignments)
  gr <- GRanges(chr, 
                IRanges(min(pos) - searchWidth, max(pos) + searchWidth), 
                strand="*"
  )

  sam <- if(subsample) {
    clippedBam(bamfile, gr, mapqFilter = mapqFilter, yieldSize = maxReads)
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
    mcols(sam)$qname  <- NULL
    
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
    mcols(sam)$qname  <- NULL

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
                            noOverClip = T,
                            searchWidth = 10L,
                            mapqFilter = 0L) {
  # to do:  subsample = T, maxReads = 1e5
  require(GenomicAlignments)
  
  gr <- GRanges(chr, 
                IRanges(min(pos) - searchWidth, max(pos) + searchWidth), 
                strand="*"
  )
  sam <- readGAlignments(bamfile,
                         param = ScanBamParam(
                           which = gr, 
                           what = c("qname", "flag", "mpos"),
                           mapqFilter = mapqFilter
                         )
  )
  sam %<>% .[isposclipped(sam, ori, pos, searchWidth = searchWidth)]
  matepos <- mcols(sam)$mpos[bamFlagTest(mcols(sam)$flag, "isProperPair")]
  if(noOverClip) {
    matepos <- if(ori == "f") {
      matepos[matepos >= pos]
    } else {
      matepos[matepos <= pos]
    }
  }
  sam %<>% .[order(njunc(.), decreasing = T)]
  sam %<>% .[!duplicated(mcols(.)$qname)]
  if(length(matepos) == 0) {
    return(list(
      pairs = GAlignmentPairs(sam[0], sam[0]),
      nopair = sam
    ))
  }
  mgr <- GRanges(chr, 
                 IRanges(min(matepos) - searchWidth, max(matepos) + searchWidth), 
                 strand="*"
  )
  mate <- readGAlignments(bamfile,
                          param = ScanBamParam(
                            which = mgr, 
                            what = c("qname", "flag"),
                            mapqFilter = mapqFilter
                          )
  ) %>% subset(qname %in% mcols(sam)$qname)
  mate %<>% .[order(njunc(.), decreasing = T)]
  
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
    paste0(stringr::str_dup(".", nchar(seq) - nchar(sseq)), sseq)
  } else {
    paste0(sseq, stringr::str_dup(".", nchar(seq) - nchar(sseq)))
  }
  seq %<>% DNAStringSet
  matchscore <- nucleotideSubstitutionMatrix(match = 0, mismatch = 1, baseOnly = TRUE) %>%
    cbind(N = 0, . = 0) %>% 
    rbind(N = 0, . = 0)
  differ <- sapply(gseq, function(x) {
    stringDist(c(seq, x), method = "substitutionMatrix", substitutionMatrix = matchscore, gapOpening = 100)
  }) / nchar(sseq)
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
                        min_length = min_length + 1) 
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
                                   searchWidth = 10L, 
                                   mapqFilter = 1L,
                                   shift_range = 0,
                                   mismatch_cutoff = 0.1, 
                                   cliplength_cutoff = 4,
                                   maxReads = 1e5,
                                   threads = getOption("mc.cores", detectCores())) {
  library(BiocParallel)
  require(GenomicAlignments)
  require(data.table)
  require(GenomicFiles)
  
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
    nonspecificTE = NA_real_
  )
  
  n <- nrow(ctea)
  if(n == 0) {
    return( lapply(NAresult, `[`, 0) )
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
    list(
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
      gap = suppressWarnings(median(meta$gap[isMatch & !shortClip])),
      secondary = mean(isSecondary[isMatch]),
      editDistance = mean(meta$NM[isMatch]),
      nonspecificTE = mean(isTEread[isMatch])
    )
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

# unique.rtea <- function(rtea, 
#                         overhang_cutoff = 10L, 
#                         maxgap = 10L, 
#                         alignscore_cutoff = 20,
#                         verbose = getOption("verbose", F)) {
#   # To do: can be parallelized over same groups
#   library(GenomicRanges)
#   library(Biostrings)
#   ungappos <- ungapPos.rtea(rtea, overhang_cutoff = overhang_cutoff)
#   rgr <- rtea[, GRanges(chr, IRanges(ungappos, ungappos), ifelse(ori == "f", "+", "-"))]
#   ovlap <- findOverlaps(rgr, rgr, maxgap = maxgap)
#   ovlap %<>% .[queryHits(.) < subjectHits(.)]
# 
#   msg("Finding duplicates...")
#   for(i in seq_along(ovlap)) {
#     if(verbose) {
#       cat(sprintf("\r%0.2f%%          ", i/length(ovlap)*100))
#     }
#     
#     queryIdx <- queryHits(ovlap[i])
#     subjectIdx <- subjectHits(ovlap[i])
#     
#     if(rtea[queryIdx, is.na(TEscore)] |
#        rtea[subjectIdx, is.na(TEscore)] |
#        rtea[queryIdx, class] != rtea[subjectIdx, class]) {
#       next
#     }
#     
#     if(rtea[c(queryIdx, subjectIdx), all(grepl("^A+$", seq)) | all(grepl("^T+$", seq))]) {
#       queryselect <- rtea[queryIdx, nchar(seq)] > rtea[subjectIdx, nchar(seq)]
#       if(rtea[queryIdx, nchar(seq)] == rtea[subjectIdx, nchar(seq)]) {
#         queryselect <- rtea[queryIdx, uniqueCnt] >= rtea[subjectIdx, uniqueCnt]
#       }
#     } else {
#       align <- pairwiseAlignment(rtea[queryIdx, seq], rtea[subjectIdx, seq], type = "local")
#       if(score(align) < alignscore_cutoff) {
#         next
#       }
#       queryselect <- rtea[queryIdx, TEscore] > rtea[subjectIdx, TEscore]
#       if(rtea[queryIdx, TEscore] == rtea[subjectIdx, TEscore]) {
#         queryselect <- rtea[queryIdx, uniqueCnt] >= rtea[subjectIdx, uniqueCnt]
#       }
#     }
# 
#     if(queryselect) {
#       rtea[subjectIdx, ] <- rtea[queryIdx, ]
#     } else {
#       rtea[queryIdx, ] <- rtea[subjectIdx, ]
#     }
#   }
#   
#   unique(rtea)
# 
# }

unique.rtea <- function(rtea, ...) {
  dup <- duplicate.rtea(rtea, ...)
  rtea[is.na(dup) | dup == seq_along(dup)]
}

duplicate.rtea <- function(rtea, 
                        overhang_cutoff = 10L, 
                        maxgap = 10L,
                        threads = 1) {
  library(GenomicRanges)
  library(Biostrings)
  opt <- options(mc.cores = threads)
  on.exit(options(opt))
  
  ungappos <- ungapPos.rtea(rtea, overhang_cutoff = overhang_cutoff)
  rgr <- rtea[, GRanges(chr, IRanges(pos, pos), ifelse(ori == "f", "+", "-"))]
  ugrgr <- rtea[, GRanges(chr, IRanges(ungappos, ungappos), ifelse(ori == "f", "+", "-"))]
  ovlap <- findOverlaps(rgr, rgr, maxgap = 10L, select = "all")
  ugovlap <- findOverlaps(ugrgr, ugrgr, maxgap = 10L, select = "all")
  ovlist <- Map(union, as.list(ovlap), as.list(ugovlap))

  lowgrp <-mclapply(ovlist, function(x) min(unlist(ovlist[x]))) %>% unlist
  grp <- mclapply(ovlist, function(x) min(lowgrp[x])) %>% unlist
  while(!identical(lowgrp, grp)) {
    lowgrp <- grp
    grp <- mclapply(ovlist, function(x) min(lowgrp[x])) %>% unlist
  }
  dupgrp <- which(table(grp) > 1) %>% names %>% as.integer
  
  bestone <- mclapply(dupgrp, function(g) {
    idx <- which(grp == g)
    highscore <- rtea[idx, TEscore] %>% {idx[. == max(., na.rm = T)]} %>% na.omit
    if(length(highscore) == 0) {
      highscore <- idx
    }
    highcnt <- rtea[highscore, uniqueCnt] %>% {highscore[. == max(.)]}
    prox <- if(rtea[highcnt, ori][1] == "f") {
      rtea[highcnt, pos] %>% {highcnt[. == max(.)]}
    } else {
      rtea[highcnt, pos] %>% {highcnt[. == min(.)]}
    }
    prox[1]
  }) %>% unlist
  bestone[match(grp, dupgrp)]
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

annotateScallop.ctea <- function(rtea, scallopfile, gencodefile, threads = getOption("mc.cores", 2)) {
  require(rtracklayer)
  msg("Loading gencode...")
  gencode <- import(gencodefile, "gtf")
  exons <- subset(gencode, type == "exon")
  rm(gencode)
  seqlevels(exons) %<>% sub("chr", "", .)
  exons %<>% .[order(.$transcript_id)]
  msg("Annotating scallop transcript...")
  
  scallop <- import(scallopfile, "gtf")
  seqlevels(scallop) %<>% sub("chr", "", .)
  
  n <- nrow(rtea)
  lsclanno <- mclapply(
    seq_len(n), 
    mc.cores = threads,
    function(idx) {
      cat(sprintf("\r%0.2f%%              ", idx/n*100))
      
      if(is.na(rtea[idx, scallop_id])) {
        return(list(scallop_type = NA))
      }
      scallop_id <- strsplit(rtea[idx, scallop_id], ",")[[1]]
      clippos <- ungapPos.rtea(rtea[idx])
      ori <- rtea[idx, ori]
      hardpos <- rtea[idx, ifelse(ori == "f", hardstart, hardend)]
      tepos <- ifelse(is.na(hardpos), clippos, hardpos)
      
      whichTranscript <- function(scl, overlapprop_cutoff = 0.7) {
        scl %<>% subset(type == "exon")
        ovtid <- subsetByOverlaps(exons, scl)$transcript_id %>% unique
        overlapsize <- sapply(ovtid, 
                              function(x) subset(exons, transcript_id == x) %>% 
                                intersect(scl) %>% 
                                width %>% 
                                sum
        ) 
        bestid <- overlapsize %>%
          .[. > sum(width(scl)) * overlapprop_cutoff] %>%
          {names(.)[which.max(.)]}
        subset(exons, transcript_id == bestid)
      }
      trpt <- whichTranscript(subset(scallop, transcript_id %in% scallop_id))
      
      if(length(trpt) == 0) {
        return(list(scallop_type = "novel transcript"))
      } else {
        strand <- as.character(strand(trpt[1]))  
        if(tepos < min(start(trpt))) {
          if(ori == "r") {
            return(list(scallop_type = NA))
          }
          if(strand == "+") {
            scallop_type <- "alternative TSS"
          } else {
            scallop_type <- "read-through"  
          }
        } else if(tepos > max(end(trpt))) {
          if(ori == "f") {
            return(list(scallop_type = NA))
          }
          if(strand == "+") {
            scallop_type <- "read-through"
          } else {
            scallop_type <- "alternative TSS"  
          }
        } else {
          scallop_type <- "exonization"
        }
        scallop_gene_id <- trpt$gene_id[1]
        scallop_gene_name <- trpt$gene_name[1]
        scallop_transcript_id <- trpt$transcript_id[1]
        scallop_transcript_name <- trpt$transcript_name[1]
        scallop_gene_type <- trpt$gene_type[1]
      }
      
      list(scallop_type = scallop_type, 
           scallop_gene_id = scallop_gene_id, 
           scallop_gene_name = scallop_gene_name,
           scallop_transcript_id = scallop_transcript_id, 
           scallop_transcript_name = scallop_transcript_name, 
           scallop_gene_type = scallop_gene_type
      )
    })
  
  data.table(rtea, rbindlist(lsclanno, fill = T))
}

fusiontype <- function(rtea, tx_id, edbpkg = "EnsDb.Hsapiens.v86") {
  stopifnot(nrow(rtea) == length(tx_id))
  
  require(edbpkg, character.only = T)
  edb <- get(edbpkg)
  
  tx <- if(sum(!is.na(tx_id)) > 0) {
    transcripts(edb,
                columns = c("seq_strand", "tx_seq_start", "tx_seq_end", "tx_biotype"),
                filter = TxIdFilter(na.omit(tx_id)),
                return.type = "GRanges") %>% data.frame %>% data.table
    # return.type = "data.frame") %>% data.table
  } else {
    data.table(tx_id = NA, start = NA, end = NA, strand = NA, tx_biotype = NA)
  }
  dt <- data.table(ugpos = ungapPos.rtea(rtea),
                   rtea[, .(ori)], 
                   tx[tx_id, on = "tx_id"]
  )
  dt[is.na(tx_id), fusion_type := "novel transcript"]
  dt[(ori == "f" & ugpos > end) | (ori == "r" & ugpos < start), fusion_type := "impossible"]
  dt[between(ugpos, start, end, NAbounds = NA), fusion_type := "exonic/exonization"]
  dt[(ori == "f" & ugpos < start), fusion_type := fifelse(strand == "+", "alternative TSS", "read-through")]
  dt[(ori == "r" & ugpos > end), fusion_type := fifelse(strand == "+", "read-through", "alternative TSS")]
  dt[, .(tx_id, tx_biotype, fusion_type)]
}

fusiontypeByCigar <- function(rtea, bamfile, 
                              edbpkg = "EnsDb.Hsapiens.v86", 
                              threads = getOption("mc.cores", detectCores())
) {
  require(BiocParallel)
  require(GenomicAlignments)
  require(ensembldb)
  
  if(paste0("package:", edbpkg) %in% search()) {
    detach(paste0("package:", edbpkg), unload = T, character.only = T)
  }
  trptMatch <- function(i, edbpkg, bamfile, chr, pos, ori) {
    cat(sprintf("\r%0.2f%%              ", i/n*100))
    
    library(edbpkg, character.only = T)
    edb <- get(edbpkg)
    
    splvec <- function(x, FUN, maxN, AGGREGATE = c, ...) {
      spl <- lapply(seq(1, length(x), maxN), function(i1) i1:min(i1 + maxN - 1, length(x)))
      lst <- lapply(spl, function(i) FUN(x[i], ...))
      do.call(AGGREGATE, lst)
    }
    
    # sam <- getClippedReads(bamfile, chr, pos, ori)
    sam <- getClippedPairs(bamfile, chr, pos, ori)
    maxlen <- 900
    
    gr <- granges(c(first(sam$pairs), second(sam$pairs), sam$nopair)) %>% unstrand
    exons <- if(length(gr) <= maxlen) {
      exons(edb, columns = "tx_id", filter = GRangesFilter(unique(gr)))
    } else {
      exons <- splvec(reduce(gr), 
                      function(x) exons(edb, columns = "tx_id", filter = GRangesFilter(x)),
                      maxN = maxlen
      )
      exons[!duplicated(data.frame(exons))]
    }
    
    ovlexon <- data.table(
      tx_id = exons$tx_id,
      cnt = findOverlaps(narrow(gr, 5, -5), exons, type = "within") %>% countSubjectHits
    ) %>% .[, .(tx_support_exon = sum(cnt)), by = tx_id]
    
    gap <- junctions(c(first(sam$pairs), second(sam$pairs), sam$nopair)) %>% unlist %>% unstrand
    
    # if(length(gap) == 0) {
    #   return(data.table(tx_id = NA_character_, tx_support_cnt = NA_integer_, numgap = length(gap)))
    # }
    gaplen <- length(unique(gap))
    if(gaplen > 0) {
      introns <- if(gaplen <= maxlen) {
        intronsByTranscript(edb, filter = GRangesFilter(unique(gap))) %>% unlist  
      } else {
        introns <- splvec(unique(gap), 
                          function(x) intronsByTranscript(edb, filter = GRangesFilter(x)),
                          maxN = maxlen
        ) %>% unlist
        # spl <- lapply(seq(1, gaplen, maxlen), function(x) x:min(x + maxlen - 1, gaplen))
        # introns <- lapply(spl, function(i) intronsByTranscript(edb, filter = GRangesFilter(unique(gap)[i]))) %>%
        #   do.call(c, .) %>% 
        #   unlist
        introns[!duplicated(data.frame(introns))]
      }
      ovlintron <- data.table(
        tx_id = names(introns),
        cnt = countOverlaps(introns, gap, type = "equal")
      ) %>% .[, .(tx_support_intron = sum(cnt)), by = tx_id]
      ovlcnt <- merge(ovlexon, ovlintron, all = T, by = "tx_id")
      ovlcnt[is.na(tx_support_exon), tx_support_exon := 0]
      ovlcnt[is.na(tx_support_intron), tx_support_intron := 0]
    } else {
      ovlcnt <- ovlexon
      ovlcnt[, tx_support_intron := 0]
    }
    
    ovlcnt %<>% .[tx_support_exon > 0 | tx_support_intron > 0]
    data.table(ovlcnt[order(-tx_support_intron, -tx_support_exon, tx_id)][1],
               numgap = length(gap))
  }
  
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
  
  ft <- fusiontype(rtea, trpt$tx_id, edbpkg = edbpkg)
  stopifnot(identical(trpt$tx_id, ft$tx_id))
  data.table(rtea, ft, trpt[, .(tx_support_exon, tx_support_intron, numgap)])
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
  ovlap %<>% .[rtea[ori == "r"][queryHits(.), ungappos] < rtea[ori == "f"][subjectHits(.), ungappos]] + maxTSD
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
  data.table(ctea, anno[, .(gene_id, gene_name, transcript_id, transcript_name, type, type_number, strand, upstream, downstream)])
  
}

localHardClip <- function(rtea,
                          mateDistMax = 100000,
                          overhangmin = 5,
                          score_cutoff = 10,
                          threads = getOption("mc.cores", detectCores())) {
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
  idxfa <- idxfa[!is.na(mfa)]
  fseq <- rtea[!is.na(mfa), DNAStringSet(seq)]
  rseq <- reverseComplement(fseq)
  ori <- rtea[!is.na(mfa), ori]
  msg("Mapping + strand to TE references...")
  n <- sum(!is.na(mfa))
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
  idx[!is.na(mfa)] <- seq_len(n)
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