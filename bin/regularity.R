#!/usr/bin/env Rscript

# Regularity - based on PSD -----------------------------------------------

# Skript by Simon Holzinger mod. by Leo Schmutterer to fit the pipeline

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(plyranges)
suppressMessages(library(zoo))
library(furrr)

args <- commandArgs(TRUE)
# Functions ---------------------------------------------------------------

# Calculate spectral power density and return the spectra
getspec <- function(x) {
  out <- spec.pgram(x, spans = c(3,3), pad = 1, plot = F)
  return(out$spec)
}

# get spectral power for each window in a rolling window approach over a coverage vector
# output is a dataframe with frequency, the corresponding NRL and the window start/end
# output is filtered to exclude small NRLs which are not relevant.
get_rolling_psd_df <- function(cov,
                               width = 513,
                               stepSize = 100,
                               lower_NRL_filter = 50) {
  start = 1
  end <- length(cov)
  
  # get power spectrum for each (rolling) window 
  specs <- (rollapply(cov[start:end],
                      width = width,
                      getspec,
                      by = stepSize))
  
  # Get frequencies 
  pgram <- spec.pgram(cov[1:width],
                      spans = c(3,3),
                      pad = 1,
                      plot = F)
  freqs <- pgram$freq
  
  colnames(specs) <- freqs
  
  # Create Dataframe
  spec_df <- specs %>% 
    as_tibble() %>% 
    mutate(start = seq(from = 1,
                       to = length(cov) - width,
                       by = stepSize),
           end = seq(from = width - 1,
                     to = length(cov),
                     by = stepSize)) %>% 
    pivot_longer(!c("start","end")) %>% 
    mutate(freq = as.numeric(name),
           NRL = 1/as.numeric(name),
           start = start + (width - 1 - stepSize) / 2,
           end = end - (width - 1 - stepSize) / 2) %>% 
    select(-name) %>% 
    filter(NRL > lower_NRL_filter) # important to lower df-size 
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
    tools::file_path_sans_ext() %>%
    str_split(pattern = "_",simplify = T) %>%
    .[1,4] %>%
    str_split(pattern = "\\.",simplify = T) %>% 
    .[1,1]
  
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
  
  return(psd_gr %>% as_tibble())
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
      #trim() %>%
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
      # trim() %>% 
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


## For parallelisation

plan(multisession, workers = 4)

# Get PSD for each timepoint
inFiles %>%
  setNames(fileNames) %>% 
  map(function(x) {
    
    print(x)
    
    # Get coverage vector
    cov_vec <- read_bigwig(x) %>% 
      coverage(weight = .$score) %>%
      lapply(as.numeric)
    
    # Get PSD over a rolling window
    y <- cov_vec %>% 
      future_map(get_rolling_psd_df,
                 width = 513,
                 stepSize = 100) %>% 
      bind_rows(.id = "seqnames")
    
    print("done...saving")
    
    # Save Intermediate Files
    name <- x %>% tools::file_path_sans_ext() %>% basename()
    saveRDS(object = y, file = paste0("PSD_RollingWindow_",name ,".rds"))
  })

# Extract track nearest to observed/expected average NRL  -----------------
avgNRL <- 180

psd_files <- list.files(pattern = "*.rds", full.names = T)

# get spectral power corresponding to the average expected NRL for each sample 
psd_df <- psd_files %>%
  map_df(get_psd_avgNRL,
         avgNRL = avgNRL,
         outPath = getwd(),
         .id = "sample") 

# calculate variance and mean from 10xLog (like in dB for sound)  --------
psd_log10_gr <- psd_df %>% 
  group_by(seqnames, start, end) %>% 
  mutate(value_log = log(value + 1) * 10) %>% 
  summarise(psd_var = var(value_log),
            psd_mean = mean(value),
            psd_sd = sd(value_log)) %>% 
  ungroup() 

# Output of variance track
var_out <- psd_log10_gr %>% 
  mutate(score = psd_var) %>% 
  filter(!is.na(score)) %>%
  as_granges()

seqlengths(var_out) <- chrSizes$size
write_bigwig(var_out, file = paste0("10xlogPSD_var_RollingWindow.bw"))

# Output of mean track
mean_out <- psd_log10_gr %>% 
  mutate(score = psd_mean) %>% 
  filter(!is.na(score)) %>%
  as_granges()

seqlengths(mean_out) <- chrSizes$size
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

cutOff<-max(which(loess.f1<3))+1
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

rankedDat %>%mutate(name = nucID)%>%
  dplyr::select(-c(width))%>%
  filter(var_rank > cutOffval) %>% 
  as_granges()%>%
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
