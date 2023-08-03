#!/usr/bin/env Rscript


#Script by Leo Schmutterer based on the "swissknife" package
#Script gets a bam file and calculates NRLs within Region defined in a bed file
args <- commandArgs(TRUE)

library(plyranges)
library(swissknife)
library(BiocParallel)
library(ggplot2)
library(plyr)

#adjusted calcPhasogramm from 'swissknife' to only load regions
calcPhasogram_regions <- function(fname, regions=NULL, rmdup=TRUE, dmax=3000L) {
  
  cnt <- numeric(dmax)
  names(cnt) <- as.character(seq.int(dmax))
  
  # get chromosomes from bam header
  bh <- Rsamtools::scanBamHeader(fname[1])[[1]]$targets
  chrs <- names(bh)
  chrsReg <- regions
  if (!is.null(regions)) {
    stopifnot(inherits(regions, "GRanges"))
    regL <- split(regions, GenomicRanges::seqnames(regions))
    chrs <- intersect(chrs, names(regL))
    regL <- regL[chrs]
  }
  
  # for each chromosome, ...
  for (i in seq_along(chrs)) {
    posPL <- posML <- list()
    for (j in seq_along(fname)) { # ... and each bam file ...
      # ... read alignment positions
      posP <- Rsamtools::scanBam(fname[j],
                                 param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE),
                                                                 what = c("pos"), which = chrsReg))[[1]]$pos
      alnM <- Rsamtools::scanBam(fname[j],
                                 param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE),
                                                                 what = c("pos","cigar"), which = chrsReg))[[1]]
      posM <- alnM$pos + GenomicAlignments::cigarWidthAlongReferenceSpace(alnM$cigar, flag = NULL, N.regions.removed = FALSE) - 1L
      posPL <- c(posPL, list(posP))
      posML <- c(posML, list(posM))
    }
    # ... combine
    posP <- do.call(c, posPL)
    posM <- do.call(c, posML)
    # ... make unique (if rmdup)
    if (rmdup) {
      posP <- unique(posP)
      posM <- unique(posM)
    }
    # ... sort
    posP <- sort(posP)
    posM <- sort(posM)
    # ... filter by selected regions
    if (!is.null(regions)) {
      posP <- posP[ IRanges::overlapsAny(IRanges::IRanges(posP, width = 1L), GenomicRanges::ranges(regL[[i]])) ]
      posM <- posM[ IRanges::overlapsAny(IRanges::IRanges(posM, width = 1L), GenomicRanges::ranges(regL[[i]])) ]
    }
    # ... phasogram (alignments on same strand)
    cnt <- calcAndCountDist(posP, posP, cnt) # plus strand
    cnt <- calcAndCountDist(posM, posM, cnt) # minus strand
  }
  
  return(cnt)
}


#set up functions for bplapply
getphasograms <- function(granges_obj) {
  result <- calcPhasogram_regions(inFile, granges_obj, )
  return(result)}

getNRL <- function(phasogram){
  estimateNRL(phasogram,
              usePeaks = 1:num_peaks_used,
              span2 = smooth_span)[1:2]
}


#setup parallel computing - ca. 4.GB memory per cpu required
num_cores <- 4
register(MulticoreParam(workers = num_cores))

#set params
smooth_span <- 0.3
num_peaks_used <- as.integer(tail(args, n = 1 ))

#read in bam
inFile <- grep(".bam$",args , value = T)

Regions <- read_bed(grep("*.bed",args,value = T))
fileNames <- args[1]
#remove brakets from names
fileNames <- gsub("[][, ]","" ,fileNames)
#split Regions into list for lapply
range_list <- split(Regions, Regions@elementMetadata@listData$name)

#calculate phasogram in regions
phasograms <- bplapply(range_list, getphasograms)

eNRL <- bplapply(phasograms, estimateNRL)

nrls <- sapply(eNRL, function(lst) lst[[1]])
#set score to nrls - NA is set to 0
Regions@elementMetadata$score <- round(nrls[as.character(Regions@elementMetadata$name)],digits = 0)
Regions@elementMetadata$score[is.na(Regions@elementMetadata$score)] <- 0

condition <- rep(fileNames, length(nrls))

nrls_values <- data.frame(nrls, condition)


write.csv(nrls_values, paste0("nrl_values_",fileNames,".csv"), quote = F, row.names = F)

write_bed(Regions, paste0("NRLs_perRegion", fileNames,".bed"))

# Create the boxplot using ggplot2
ylim1<- boxplot.stats(nrls_values$nrls)$stats[c(1,5)]

meds <- ddply(nrls_values, .(condition), summarise, med = round(median(nrls, na.rm = T),digits=2))

pdf(paste0("Boxplot_of_NRL_in_regions.pdf"))
ggplot(nrls_values, aes(x = condition, y = nrls)) +
  geom_boxplot()+
  coord_cartesian(ylim = ylim1*1.05)+
  labs(title = paste0("NRLs ", fileNames))+ xlab(condition)+ ylab("NRLs")+theme(axis.text.x = element_blank())+
  geom_text(data = meds, aes(x = condition, y = med, label = med), 
            size = 3, vjust = -0.5)
dev.off()




