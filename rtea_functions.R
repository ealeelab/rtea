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
sacct <- function(start_date){
  Sopt <- if(missing(start_date)) {
    NULL
  } else {
    paste("-S", start_date)
  }
  saccmd <- paste("sacct -P -u $USER --format=JobName%30,Partition,state,ExitCode,NCPUS,Elapsed,End,MaxRSS --units=G",
                  Sopt,
                  "| grep -v '\\---'"
  )
  sac <- fread(saccmd,
               colClasses = c(MaxRSS = "character"),
               fill = T
  )
  batchidx <- which(sac$JobName == "batch")
  sac[batchidx - 1, ExitCode := sac$ExitCode[batchidx]]
  sac[batchidx - 1, MaxRSS := sac$MaxRSS[batchidx]]
  sac <- sac[-batchidx]
  sac[, Elapsed := as.difftime(Elapsed, units = "hours")]
  ystday <- Sys.Date() - 1
  plocale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  sac[grep("Ystday", End), End := sub("Ystday", format(ystday, "%d %b"), End)]
  Sys.setlocale("LC_TIME", plocale)
  endformat <- rep("%d %b %H:%M", nrow(sac))
  endformat[grep(" ", sac$End, invert = T)] <- "%H:%M:%S"
  suppressWarnings(sac[, End := strptime(End, endformat)])
  sac
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

readctea <- function(cteafile, count_cutoff = 3, confidence_cutoff = 2, threads = getDTthreads()) {
  require(Biostrings)
  require(parallel)
  setDTthreads(threads)
  writeLines(paste("Reading", cteafile))
  # ctea <- fread(cteafile)
  awkcmd <- paste("awk -F '\\t' '$4 >=", count_cutoff, 
                  "&& $6 >=", confidence_cutoff, 
                  "'",
                  # "&& $5 != \"PolyA\"'", 
                  cteafile)
  ctea <- fread(cmd = awkcmd, sep = "\t")
  setnames(ctea, 
           c("chr", "pos", "ori", "cnt", "family", "confidence",
             "moreFamily", "seq", "morePos", "moreSeq"))
  # ctea <- ctea[cnt >= count_cutoff & confidence >= confidence_cutoff, ]
  revseq <- ctea[["seq"]] %>% 
    DNAStringSet %>% 
    reverseComplement %>% 
    as("character")
  moreseq <- strsplit(ctea[["moreSeq"]], ",")
  morepos <- strsplit(ctea[["morePos"]], ",")
  seqmatch <- mcmapply(match, ctea[["seq"]], moreseq, mc.cores = threads)
  revseqmatch <- mcmapply(match, revseq, moreseq, mc.cores = threads)
  seqmatch[is.na(seqmatch)] <-revseqmatch[is.na(seqmatch)]
  ctea[, pos := mcmapply(`[`, morepos, seqmatch, mc.cores = threads) %>% as.integer]
  ctea[, seq := mcmapply(`[`, moreseq, seqmatch, mc.cores = threads)]
  ctea[, class := NA_character_]
  ctea[family %in% c("L1HS", "LINE1"), class := "L1"]
  ctea[toupper(substr(family, 1, 3)) == "ALU", class := "Alu"]
  ctea[substr(family, 1, 3) == "SVA", class := "SVA"]
  ctea[substr(family, 1, 4) == "HERV", class := "HERV"]
  # ctea[, class := structure(c(rep("L1", 2),
  #                             rep("Alu", 11),
  #                             rep("SVA", 7),
  #                             rep("HERV", 8),
  #                             rep("simple_repeat", 3)),
  #                           names = c("L1HS", "LINE1", 
  #                                     "ALU", "AluJb", "AluJo", "AluSc", "AluSg", "AluSp", "AluSq", "AluSx", "AluSz", "AluY", "AluYd8", 
  #                                     "SVA", "SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F",
  #                                     "HERVE", "HERVH48I", "HERVK", "HERVK11DI", "HERVK11I", "HERVK13I", "HERVK3I", "HERVK9I", 
  #                                     "(A)n", "(GAAA)n", "(TTTC)n"
  #                           ))[family]]
  ctea[, .(chr, pos, ori, cnt, class, family, moreFamily, confidence, seq)]
}

readxtea <- function(xteafile, total_count_cutoff = 5, threads = getDTthreads()) {
  setDTthreads(threads)
  writeLines(paste("Reading", xteafile))
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
  matchdata <- matchdata[, lapply(.SD, function(x) paste(unique(x), collapse = ",")), by = idx]
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

filterTEsite.ctea <- function(ctea, clip_width = 5,
                              filterFamily = c("Alu", "Simple_repeat", "Low_complexity")
                     ) {
  rpt <- ctea[, repeatLocation(chr, pos - clip_width, pos + clip_width)]
  ctea[, c("posRepFamily", "posRep") := rpt]
  ctea <- ctea[strsplit(posRepFamily, ",") %>%
         sapply(function(x) ! any(x %in% filterFamily)), ]
  ctea[!grepl("SVA", family) | !grepl("SVA", posRep)]
}

filterPolyTEsite.ctea <- function(ctea, clip_width = 5) {
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

filterPolyA.ctea <- function(ctea, polyAoutfile = NULL, Alength = 10, percentA = 0.8) {
  library(stringr)
  require(GenomicRanges)
  
  normalPolyA <- rep(NA, nrow(ctea))
  normalPolyA[ctea[, ori == "r"]] <- ctea[ori == "r", 
                                          grepl(paste0("A{", Alength, "}"), seq) | 
                                            str_count(seq, "A") / nchar(seq) > percentA
                                          ]
  normalPolyA[ctea[, ori == "f"]] <- ctea[ori == "f",
                                          grepl(paste0("T{", Alength, "}"), seq) | 
                                            str_count(seq, "T") / nchar(seq) > percentA
                                          ]

  ctea[normalPolyA == F, ]
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

filterNoClip.ctea <- function(ctea, similarity_cutoff = 75, threads = getDTthreads()) {
  require(bsgenomePKG,  character.only = T)
  genome <- get(bsgenomePKG)
  ctea[, start := pos]
  ctea[, end := pos]
  if(nrow(ctea) == 0) {
    writeLines("\nNo insertion found.")
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
  writeLines("Comparing clipped sequences with reference sequences...")
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
  writeLines("\nComparing done.")
  ctea[, similar := unlist(similar)]
  ctea[similar <= similarity_cutoff, ]
}

getClippedReads <- function(bamfile, chr, pos, ori = c("f", "r"), 
                            searchWidth = 10L,
                            mapqFilter = 10L) {
  require(GenomicAlignments)
  gr <- GRanges(chr, 
                IRanges(min(pos) - searchWidth, max(pos) + searchWidth), 
                strand="*"
  )
  sam <- readGAlignments(
    bamfile, 
    param = ScanBamParam(
      which = gr, 
      what = c("qname", "seq", "qual", "flag", "mpos"),
      tag = "NM",
      mapqFilter = mapqFilter
    )
  ) 
  
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
    split <- grepl("^[0-9]+S", cigar(sam))
    othersplit <- grepl("S$", cigar(sam))
    sam <- sam[order(secondary, split, othersplit)]
    isFirst <- bamFlagTest(mcols(sam)$flag, "isFirstMateRead")
    sam <- sam[!duplicated(paste(mcols(sam)$qname, isFirst)) & grepl("^[0-9]+S", cigar(sam))]
    
    mcols(sam)$slen <- sub("S.*", "", cigar(sam)) %>% as.integer
    mcols(sam)$shift <- start(sam) - pos - 1
    if(length(sam) == 0) {
      mcols(sam) <- assignZero(mcols(sam))
      return(sam)
    }
    sam <- sam[abs(mcols(sam)$shift) < searchWidth & 
                 mcols(sam)$slen - mcols(sam)$shift > 0]
    mcols(sam)$sseq <- subseq(mcols(sam)$seq, 1, mcols(sam)$slen - mcols(sam)$shift)
    mcols(sam)$squal <- subseq(mcols(sam)$qual, 1, mcols(sam)$slen - mcols(sam)$shift) %>%
      as("IntegerList") %>%
      mean
    mcols(sam)$overhang <- sub("^[0-9]+S([0-9]+)M[0-9]+N.*", "\\1", cigar(sam)) %>%
      as.integer
    mcols(sam)$gap <- sub("^[^N]+?([0-9]+)N.*", "\\1", cigar(sam)) %>%
      as.integer
    # sam <- sam[order(mcols(sam)$squal, decreasing=T)]
    # sam <- sam[!duplicated(paste(end(sam), mcols(sam)$mpos))]
    # sam <- sam[!duplicated(paste(mcols(sam)$sseq, strand(sam)))]
    
  } else if(ori == "r") {
    # deduplication
    sam <- sam[end(sam) < pos + qwidth(sam)]
    secondary <- bamFlagTest(mcols(sam)$flag, "isSecondaryAlignment")
    split <- grepl("S$", cigar(sam))
    othersplit <- grepl("^[0-9]+S", cigar(sam))
    sam <- sam[order(secondary, split, othersplit)]
    isFirst <- bamFlagTest(mcols(sam)$flag, "isFirstMateRead")
    sam <- sam[!duplicated(paste(mcols(sam)$qname, isFirst)) & grepl("S$", cigar(sam))]
    
    # sam <- sam[grep("^[0-9]+S", cigar(sam), invert = T)]
    mcols(sam)$slen <- 
      explodeCigarOpLengths(cigar(sam), ops="S") %>%
      sapply(function(x) x[length(x)])
    mcols(sam)$shift <- end(sam) - pos + 1
    if(length(sam) == 0) {
      mcols(sam) <- assignZero(mcols(sam))
      return(sam)
    }
    sam <- sam[abs(mcols(sam)$shift) < searchWidth & 
                 mcols(sam)$slen + mcols(sam)$shift > 0]
    mcols(sam)$sseq <- subseq(mcols(sam)$seq, -(mcols(sam)$slen + mcols(sam)$shift))
    mcols(sam)$squal <- subseq(mcols(sam)$qual, -(mcols(sam)$slen + mcols(sam)$shift)) %>%
      as("IntegerList") %>%
      mean
    mcols(sam)$overhang <- sub(".*N([0-9]+)M[0-9]+S$", "\\1", cigar(sam)) %>%
      as.integer
    mcols(sam)$gap <- sub(".*?([0-9]+)N[^N]+$", "\\1", cigar(sam)) %>%
      as.integer
    # sam <- sam[order(mcols(sam)$squal, decreasing=T)]
    # sam <- sam[!duplicated(paste(start(sam), mcols(sam)$mpos))]
    # sam <- sam[!duplicated(paste(mcols(sam)$sseq, strand(sam)))]
  }
  
  sam
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
  seq <- DNAStringSet(seq)
  differ <- lapply(sseq, function(x) {
    subseqprox(seq, nchar(x)) %>%
      c(DNAStringSet(x)) %>%
      stringDist %>%
      `/`(nchar(x))
  }) %>% unlist
  
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
                           falseCntRatioCutoff = 3,
                           baseQualCutoff = 25,
                           wrongPosPropCutoff = 0.2,
                           bothClipPropCutoff = 0.2,
                           secondaryCutoff = 0.99,
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
                                   searchWidth = 10, 
                                   shift_range = 0,
                                   mismatch_cutoff = 0.1, 
                                   cliplength_cutoff = 4,
                                   threads = getDTthreads()) {
  library(parallel)
  require(GenomicAlignments)
  require(data.table)
  n <- nrow(ctea)
  if(n == 0) {
    return( data.table(
      ctea, 
      matchCnt = integer(0),
      trueCnt = integer(0), 
      uniqueCnt = integer(0),
      bothClip = integer(0),
      polyAcnt = integer(0),
      falseCnt = integer(0), 
      discCnt = integer(0),
      baseQual = numeric(0), 
      clustered = numeric(0),
      anyOverClip = integer(0),
      mateDist = integer(0),
      anyWrongPos = integer(0),
      overhang = integer(0),
      gap = integer(0),
      secondary = numeric(0),
      editDistance = numeric(0),
      nonspecificTE = numeric(0)
    ) )
  }
  writeLines("Counting clipped reads...")
  lcnt <- mcmapply(
    seq_len(n), 
    bamfile,
    ctea$chr, ctea$pos, ctea$ori, ctea$seq, ctea$family,
    mc.cores = threads, 
    mc.preschedule = T, 
    SIMPLIFY = F,
    FUN = function(i, bamfile, chr, pos, ori, seq, family) {
      if(n > 10) cat(sprintf("\r%0.2f%%              ", i/n*100))
      if(anyNA(c(chr, pos, ori, seq))) {
        return(list(
          matchCnt = NA_integer_,
          trueCnt = NA_integer_, 
          uniqueCnt = NA_integer_,
          bothClip = NA_integer_,
          polyAcnt = NA_integer_,
          falseCnt = NA_integer_, 
          discCnt = NA_integer_,
          baseQual = NA_real_, 
          clustered = NA_real_,
          anyOverClip = NA_integer_,
          mateDist = NA_integer_,
          anyWrongPos = NA_integer_,
          overhang = NA_integer_,
          gap = NA_integer_,
          secondary = NA_real_,
          editDistance = NA_real_,
          nonspecificTE = NA_real_
        ))
      }
      sam <- getClippedReads(bamfile, chr, pos, ori, searchWidth)
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
      isTEread <- TEalignScore(mcols(sam)$seq, family)
      list(
        matchCnt = sum(isMatch),
        trueCnt = sum(isMatch & !bothClip & !shortClip & !isPolyA), 
        uniqueCnt = sum(isMatch & !possibleDup & !shortClip),
        bothClip = sum(isMatch & bothClip),
        polyAcnt = sum(isPolyA),
        falseCnt = sum(!isMatch & !shortClip), 
        discCnt = sum(isMatch & !isProperPair & !mateUnmapped & isMateSide),
        baseQual = median(meta$squal[isMatch & !isPolyA]), 
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
  )
  writeLines("\nCounting done.")
  tryCatch(
    cntdt <- rbindlist(lcnt),
    error = function(e) {
      # nelem <- sapply(lcnt, function(x) length(unlist(x)))
      # if(any(nelem != 15)) {
      #   message("Missing element in row number: ", 
      #           toString(which(nelem != 15))
      #   )
      ismsg <- sapply(lcnt, class) == "character"
      if(any(ismsg)) {
        message(lcnt[ismsg] %>% unlist %>%  unique %>% paste(collapse = "\n"))
      }
      e
    }
  )
  
  data.table(ctea, cntdt)
}

countClippedReads.default <- function(chr, pos, ori, seq, 
                                      bamfile, 
                                      searchWidth = 10, 
                                      shift_range = 0,
                                      mismatch_cutoff = 0.1, 
                                      cliplength_cutoff = 4,
                                      threads = getDTthreads()) {
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
                                  threads = getDTthreads()) {
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

nearbyfusions <- function(ctea, matchrange = 10000, maxTSD = 10, overhangmin = 5) {
  library(GenomicRanges)
  rctea <- ctea[ori == "r"]
  fctea <- ctea[ori == "f"]
  rctea[overhang <= overhangmin, pos := pos - gap]
  fctea[overhang <= overhangmin, pos := pos + gap]
  rgr <- rctea[, GRanges(chr, IRanges(pos - maxTSD, pos + matchrange), strand = "*")]
  fgr <- fctea[, GRanges(chr, IRanges(pos, pos), strand = "*")]

  ovl <- findOverlaps(rgr, fgr)
  sameclass <- rctea[queryHits(ovl), class] == fctea[subjectHits(ovl), class] |
    xor(rctea[queryHits(ovl), isPolyA], fctea[subjectHits(ovl), isPolyA])
  nearctea <- rbind(rctea[queryHits(ovl[sameclass])],
                    fctea[subjectHits(ovl[sameclass])])
  unique(nearctea[order(as.integer(chr), pos)])
  
}

geneticLocation <- function(chr, start, end = start) {
  # require(rtracklayer)
  # writeLines(paste("Importing", gene_file))
  # annodata <- import(gene_file, "gtf")
  annodata <- readRDS(gene_file)
  writeLines("Location idenfying")
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
  annodata <- readRDS(gene_file)
  seqlevels(annodata) %<>% pastechr
  annodata$exon_number <- as.integer(annodata$exon_number)
  writeLines("Annotating TEI.")
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
                          threads = getDTthreads()) {
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
  writeLines("Hard clip align ...")
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
  writeLines("\nDone align.")
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

TEcoordinate <- function(rtea, threads = getDTthreads()) {
  require(Biostrings)
  require(parallel)
  writeLines("Reading TE fasta file...")
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
  writeLines("Mapping + strand to TE references...")
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
  writeLines("\nMapping - strand to TE references ...")
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
  writeLines("\nMapping done.")
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

readrteas <- function(rteafiles) {
  rteas <- lapply(rteafiles, fread)
  names(rteas) <- basename(rteafiles) %>% sub(".rnatea.txt", "", .)
  rteas <- rbindlist(rteas, idcol = "ID", use.names = T)
  rteas[, confident := F]
  rteas[trueCnt > 2 & falseCnt < 10 & anyOverClip == F & baseQual > 25, confident := T]
  rteas[trueCnt == 3 & falseCnt > 2, confident := F]
  cnt <- rteas[, .N, .(chr, pos, ori, gene_name, type_number)]
  merge(cnt, rteas)
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