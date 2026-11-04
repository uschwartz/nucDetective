#!/usr/bin/env Rscript

# Regularity - based on PSD -----------------------------------------------

# Skript by Simon Holzinger mod. by Leo Schmutterer to fit the pipeline

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(plyranges)
suppressMessages(library(zoo))
library(furrr)

args <- commandArgs(TRUE)
# extract average NRL from args
avgNRL <- as.numeric(args[length(args)])
# remove the last argument (avgNRL) from args
args <- args[-length(args)]

# Functions ---------------------------------------------------------------

# Calculate spectral power density and return the spectra.
# `keep` optionally restricts the returned spectrum to a subset of frequency indices
# (the periodogram FFT is computed in full regardless; only the returned vector shrinks).
getspec <- function(x, keep = NULL) {
  out <- spec.pgram(x, spans = c(3,3), pad = 1, plot = F)
  if (is.null(keep)) return(out$spec)
  return(out$spec[keep])
}

# get spectral power for each window in a rolling window approach over a coverage vector
# output is a dataframe with frequency, the corresponding NRL and the window start/end
# output is filtered to exclude small NRLs which are not relevant.
get_rolling_psd_df <- function(cov,
                               width = 513,
                               stepSize = 100,
                               lower_NRL_filter = 50,
                               avgNRL = NULL) {
  start = 1
  end <- length(cov)

  # Get frequencies (identical for every window of this width)
  pgram <- spec.pgram(cov[1:width],
                      spans = c(3,3),
                      pad = 1,
                      plot = F)
  freqs <- pgram$freq

  # Only frequencies whose NRL (= 1/freq) exceeds the filter are ever kept downstream
  # (filter(NRL > lower_NRL_filter)). Restrict the spectrum to those *before* the
  # rolling apply / pivot, so the intermediate is ~25x smaller for a large genome
  # instead of being built in full and then discarded. Rows returned are identical.
  keep <- which((1 / as.numeric(as.character(freqs))) > lower_NRL_filter)

  # When avgNRL is supplied, restrict further to the frequency/frequencies whose
  # rounded NRL is nearest avgNRL - that (and only that) is what get_psd_avgNRL
  # keeps downstream, and the frequency grid is identical for every window, so
  # selecting here is equivalent. Doing it BEFORE the rolling apply / pivot
  # collapses the per-chromosome intermediate from ~20 freqs/window to 1, which
  # is the memory pole on a large genome.
  if (!is.null(avgNRL)) {
    nrl_r <- round(1 / as.numeric(as.character(freqs))[keep])
    keep  <- keep[abs(nrl_r - avgNRL) == min(abs(nrl_r - avgNRL))]
  }
  freqs <- freqs[keep]

  # Power spectrum of every rolling window (windows: cov[s : s+width-1] for
  # s = 1, 1+stepSize, ...). Was zoo::rollapply(..., by = stepSize); over a whole
  # chromosome that piled up ~20 GB of transient spec.pgram garbage before R's GC
  # swept it - millions of calls with no breathing room. Same windows, computed
  # in blocks with a gc() between them, keeps the peak near 1 GB. Verified
  # byte-identical to rollapply (scratchpad/).
  win_starts <- seq.int(1L, length(cov) - width + 1L, by = stepSize)
  nkeep <- length(keep)
  specs <- matrix(0, nrow = length(win_starts), ncol = nkeep)
  block <- 50000L
  for (b0 in seq.int(1L, length(win_starts), by = block)) {
    ix <- b0:min(b0 + block - 1L, length(win_starts))
    chunk <- vapply(win_starts[ix],
                    function(s) getspec(cov[s:(s + width - 1L)], keep),
                    numeric(nkeep))
    specs[ix, ] <- if (nkeep == 1L) chunk else t(chunk)
    gc(FALSE)
  }
  colnames(specs) <- freqs

  # Create Dataframe.
  # start/end are taken from nrow(specs): windows run to the last full window,
  # which for some sequence lengths is one more than the old
  # seq(from = 1, to = length(cov) - width, by = stepSize) produced - that
  # mismatch threw a size error in mutate(). Identical values otherwise.
  spec_df <- specs %>%
    as_tibble() %>%
    mutate(start = seq(from = 1,
                       by = stepSize,
                       length.out = nrow(specs)),
           end = seq(from = width - 1,
                     by = stepSize,
                     length.out = nrow(specs))) %>%
    pivot_longer(!c("start","end")) %>%
    mutate(freq = as.numeric(name),
           NRL = 1/as.numeric(name),
           start = start + (width - 1 - stepSize) / 2,
           end = end - (width - 1 - stepSize) / 2) %>%
    select(-name) %>%
    filter(NRL > lower_NRL_filter) # no-op after the restriction above; kept for clarity
  return(spec_df)
}


# get spectral power corresponding to the expected average NRL. 
# Saves a bigwig file with the spectral power as score.
# returns a tibble with the scores and windows
get_psd_avgNRL <- function(x, avgNRL, outPath) {
  # read PSD data 
  psd_dat <- readRDS(x)
  # get name of current time point
  ctp <- x %>% basename() %>% 
    tools::file_path_sans_ext() 
  
  # Round and filter for NRL nearest to avg NRL
  psd_gr <- psd_dat %>%
    mutate(NRL = round(NRL),
           score = value) %>%
    filter(abs(NRL - avgNRL) == min(abs(NRL - avgNRL))) %>% 
    as_granges()
  
  # Correct for different seqlevels
  corrected_chrSizes <- chrSizes %>% filter(seqnames %in% seqnames(psd_gr))
  
  seqlengths(psd_gr) <- corrected_chrSizes$size
  write_bigwig(psd_gr, file = paste0("PSD_avgNRL_",avgNRL,"_", ctp, ".bw"))

  # only these four columns are used downstream (psd_log10_gr groups by
  # seqnames/start/end and aggregates value); returning the full ~9-column tibble
  # for every sample makes the stacked psd_df twice as large for no reason.
  return(psd_gr %>% as_tibble() %>% dplyr::select(seqnames, start, end, value))
}


## Function from https://divingintogeneticsandgenomics.rbind.io/post/compute-averages-sums-on-granges-or-equal-length-bins/
## Added Support for Rolling Bins 
averagePerBin <- function(x, binsize, mcolnames=NULL, overlap = NULL, ref_positions = NULL)
{
  if (!is(x, "GenomicRanges"))
    stop("'x' must be a GenomicRanges object")
  if (any(is.na(seqlengths(x))))
    stop("'seqlengths(x)' contains NAs")
  # If reference regions are provided we use them centered with binsize as bins
  if (!is.null(ref_positions)) {
    bins <- ref_positions %>% 
      anchor_center() %>%           # center regions
      mutate(width = binSize) %>%   # extend to binsize
      sortSeqlevels() %>%           # sort
      sort()
    seqlevels(bins) <- seqlevels(x)
    seqinfo(bins) <- seqinfo(x)
    bins <- bins %>% 
      trim() %>%
      as_tibble() %>% 
      split(f = .$seqnames) %>%
      map(~as_iranges(.)) %>%
      IRangesList(.)
  } else if (overlap >= 1) {
    tiles <- tileGenome(seqinfo(x),
                        tilewidth = overlap,
                        cut.last.tile.in.chrom = T)
    windows <- resize(tiles, binsize) # you will get a warning about trimming
    bins <- windows %>% 
      trim() %>% 
      as_tibble() %>%
      split(f = .$seqnames) %>%
      map(~as_iranges(.)) %>%
      IRangesList(.)
  } else {
    bins <- IRangesList(lapply(seqlengths(x),
                               function(seqlen)
                                 IRanges(breakInChunks(seqlen, chunksize = binsize))))
    
  }
  
  ans <- as(bins, "GRanges")
  seqinfo(ans) <- seqinfo(x)
  if (is.null(mcolnames))
    return(ans)
  averageMCol <- function(colname)
  {
    cvg <- coverage(x, weight = colname)
    views_list <- RleViewsList(
      lapply(names(cvg),
             function(seqname)
               Views(cvg[[seqname]], bins[[seqname]])))
    unlist(viewMeans(views_list), use.names = FALSE)
  }
  mcols(ans) <- DataFrame(lapply(mcols(x)[mcolnames], averageMCol))
  ans
}

# Calculate PSD -----------------------------------------------------------

## Set directories
inFiles <- grep(".bw",args , value = T)

fileNames <- args[-grep(".bw",args)]
fileNames <- fileNames[-grep(".txt", fileNames)]
fileNames <- fileNames[-grep(".bed", fileNames)]
fileNames <- gsub("[][, ]","" ,fileNames)

#optional ref_positions
if(length(which(grepl(".bed", args)) == T) == 1){
  ref_positions <- read_bed(grep(".bed", args, value = T))
  }else{ref_positions <- NULL}


chrSizes <- read.table(grep(".txt",args, value = T))
colnames(chrSizes) <- c("seqnames","size")


## PSD per chromosome, run in parallel across cores. Was plan(multicore,
## workers = 4) + future_map; the forked workers died inside the container
## ("forked worker no longer alive" - Docker's small default /dev/shm).
## plan(multisession) uses PSOCK workers (sockets, no fork, no /dev/shm), and
## each worker reads only its own chromosome from the bigwig so just one
## chromosome-sized vector crosses to it. Worker count from REGULARITY_WORKERS
## (the module sets it to task.cpus); the fallback is deliberately small so a
## bare manual run does not spawn one heavy worker per host core. spec.pgram
## uses no RNG, so the result is identical to sequential (verified in scratchpad/).
workers <- suppressWarnings(as.integer(Sys.getenv("REGULARITY_WORKERS", "")))
if (is.na(workers) || workers < 1L) workers <- min(4L, max(1L, future::availableCores()))
plan(multisession, workers = workers)

psd_one_chr <- function(chr, bwFile, chrLenAll, avgNRL) {
  chrLen <- chrLenAll[[chr]]
  gr <- rtracklayer::import.bw(
          bwFile,
          which = GenomicRanges::GRanges(chr, IRanges::IRanges(1, chrLen)))
  if (length(gr) == 0L) {
    cc <- S4Vectors::Rle(0, chrLen)
  } else {
    cc <- GenomicRanges::coverage(gr, weight = gr$score)[[chr]]
  }
  get_rolling_psd_df(as.numeric(cc), width = 513, stepSize = 100, avgNRL = avgNRL)
}

# Get PSD for each timepoint
inFiles %>%
  setNames(fileNames) %>%
  map(function(x) {

    print(x)

    chrLen <- GenomeInfoDb::seqlengths(seqinfo(rtracklayer::BigWigFile(x)))

    # data (bwFile, chrLenAll, avgNRL) is forwarded to every call and exported
    # automatically; only the free helper functions need listing in globals.
    # chunk_size = 1: one future per chromosome. Without it furrr statically
    # splits the chromosomes into one block per worker up front - and since the
    # bigwig lists the big chromosomes first, one worker gets chr1 + chr10..16
    # (~40 % of the work) and runs them serially while the rest idle. One future
    # per chromosome lets a free worker pick up the next chromosome immediately.
    y <- names(chrLen) %>%
      set_names() %>%
      future_map(
        psd_one_chr,
        bwFile = x, chrLenAll = chrLen, avgNRL = avgNRL,
        .options = furrr_options(
          seed = TRUE,
          chunk_size = 1,
          globals = c("psd_one_chr", "get_rolling_psd_df", "getspec"),
          packages = c("rtracklayer", "GenomicRanges", "IRanges",
                       "GenomeInfoDb", "S4Vectors", "zoo",
                       "tibble", "dplyr", "tidyr"))) %>%
      bind_rows(.id = "seqnames")

    print("done...saving")

    # Save Intermediate Files
    name <- x %>% tools::file_path_sans_ext() %>% basename()
    saveRDS(object = y, file = paste0("PSD_RollingWindow_", name, ".rds"))
  })

# The PSOCK workers are done with; shut them down so their ~1-2 GB each of loaded
# packages is not sitting idle through the serial aggregation below (that plus the
# grouped summarise plus averagePerBin's coverage was enough to OOM a 72 GB run).
plan(sequential)

# Extract track nearest to observed/expected average NRL  -----------------


psd_files <- list.files(pattern = "*.rds", full.names = T)

# get spectral power corresponding to the average expected NRL for each sample
psd_df <- psd_files %>%
  map_df(get_psd_avgNRL,
         avgNRL = avgNRL,
         outPath = getwd(),
         .id = "sample")

# calculate variance and mean from 10xLog (like in dB for sound)  --------
# log(value+1)*10 is evaluated inside summarise() rather than added as a grouped
# mutate() column first - same numbers, but avoids materialising a 124 M-row
# value_log column (and the grouped-mutate pass) on a human-genome run.
psd_log10_gr <- psd_df %>%
  group_by(seqnames, start, end) %>%
  summarise(psd_var  = var(log(value + 1) * 10),
            psd_mean = mean(value),
            psd_sd   = sd(log(value + 1) * 10),
            .groups  = "drop")

# Output of variance track
var_out <- psd_log10_gr %>% 
  mutate(score = psd_var) %>% 
  filter(!is.na(score)) %>%
  as_granges()

#extract correct chrom sizes
chrSizes.match<-chrSizes$size[match(names(seqlengths(var_out)),chrSizes$seqnames)]

seqlengths(var_out) <- chrSizes.match
write_bigwig(var_out, file = paste0("10xlogPSD_var_RollingWindow.bw"))

# Output of mean track
mean_out <- psd_log10_gr %>% 
  mutate(score = psd_mean) %>% 
  filter(!is.na(score)) %>%
  as_granges()

seqlengths(mean_out) <- chrSizes.match
write_bigwig(mean_out, file = paste0("10xlogPSD_mean_RollingWindow.bw"))



# Find regions with high variance in regularity ---------------------------

binSize <- 800
cutoff_slope <- 3

# Calculate average variance for overlapping bins and rank
rankedDat <- averagePerBin(var_out,
                           binsize = binSize,
                           mcolnames = "score",
                           overlap = 200, ref_positions = ref_positions) %>%
  as_tibble()

  rankedDat <- cbind(rankedDat,if (!is.null(ref_positions)){(nucID <- paste0("nuc",1:length(rankedDat$start)))}else{paste0("bin",1:length(rankedDat$start))}) %>%
  group_by(seqnames, start, end) %>% 
  ungroup() %>% 
  filter(width > 500) %>% #bins are trimmed at the ends and starts of chromosomes
  arrange(score) %>% 
  mutate(var_rank = 1:n(),
         var_rank = var_rank/max(var_rank),
         norm_var = score/max(score))  

if (!is.null(ref_positions)){
colnames(rankedDat)[7]<-"nucID"}else{colnames(rankedDat)[7]<-"bins"}


span.loess.pre<-1000/length(rankedDat$norm_var)
span.loess<-ifelse(span.loess.pre<0.05, span.loess.pre, 0.05)

lm<-loess(rankedDat$norm_var~rankedDat$var_rank,span=span.loess)
#Curve fitting
loess.line<-predict(lm)

#calculate the first derivative
loess.f1<-diff(loess.line)/diff(rankedDat$var_rank)

cutOff<-max(which(loess.f1<cutoff_slope))+1
cutOffval<- rankedDat$var_rank[cutOff]
# Plot cutoff
png(paste0("bs", binSize,"_varRegularity_selection.png"),width = 1280, height = 1280, res =300)
plot(rankedDat$var_rank, rankedDat$norm_var,
     ylab = "norm 10xLog of reference position", 
     xlab = "norm rank", bty = "n", pch = 19)
points(rankedDat$var_rank[cutOff:length(rankedDat$norm_var)],
       rankedDat$norm_var[cutOff:length(rankedDat$norm_var)],
       pch=19, col="#7F81B0")
lines(rankedDat$var_rank,loess.line, col="red")
abline(v = rankedDat$var_rank[cutOff], col = "blue", lty = 2)
dev.off()

# select the peaks with highest variance and write beds

rankedDat %>%
  mutate(name = ifelse(!is.null(ref_positions),nucID,bins))%>%
  dplyr::select(-c(width))%>%
  filter(var_rank > cutOffval) %>% 
  arrange(desc(norm_var))%>%as_granges()%>%
  write_bed(file = paste0("Irregular_Regions_bs",binSize, "_10xlog.bed"))




rankedDat %>%
  filter(var_rank > cutOffval) %>% 
  as_granges() %>%
  reduce() %>% # collapse bins into single bin if they are overlapping
  write_bed(file = paste0("Irregular_Regions_bs", binSize, "_10xlog_reduced_filt.bed"))


colnames(rankedDat)[1]<-"chr"
category <- rep("normal", length(rankedDat$var_rank))
category[which(rankedDat$var_rank> cutOffval)]<-"dynamicRegularity"
rankedDat<-cbind(rankedDat,category)
rankedDat %>% dplyr::select(-c(width, score, var_rank))%>%arrange(desc(norm_var))%>%
write.table(file="Regularity_10xlog_result_table.tsv",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = T)
