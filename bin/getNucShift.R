#!/usr/bin/env Rscript
library(rtracklayer)
library(parallel)

# Detect nucleosome shifts: for every reference nucleosome and every sample, fit a
# loess curve to the bigWig signal within the nucleosome and take the position of
# the highest local maximum; nucleosomes whose summit position varies most between
# samples are reported.
#
# Speed-up (results unchanged): each bigWig is read once for all reference
# nucleosomes instead of once per nucleosome, and the per-nucleosome loess fits run
# in parallel on SHIFT_WORKERS cores (set by the Nextflow module from task.cpus).
#
# Memory: the per-nucleosome signal is kept as a few flat vectors per sample plus
# an offset index, not as lists with one small vector per nucleosome. With ~10 M
# reference nucleosomes such lists hold tens of millions of R objects, and every
# forked worker's garbage collector writes to all of them, so each worker ended up
# with its own copy of the whole signal (copy-on-write). Flat vectors stay shared.

# load reference nucleosomes top20
nucRef<-import.bed("NucPosRef_top20.bed")


files.bw<-grep("bw",list.files(), value = T)
#get name
names.bw<-gsub("_qNorm.bw","",files.bw)

workers <- suppressWarnings(as.integer(Sys.getenv("SHIFT_WORKERS", "1")))
if (is.na(workers) || workers < 1) workers <- 1


n_nuc <- length(nucRef)

## read every bigWig once and assign its records to the reference nucleosomes.
## import.bw(which = <one range>) clips records to that range, so the records
## are clipped to their nucleosome here as well (identical input to the loess fits).
## findOverlaps() returns the hits sorted by nucleosome, so the records of
## nucleosome i are pos/score[offset[i] + 1 .. offset[i + 1]], in subject order.
signal_per_nuc <- lapply(files.bw, function(f){
    sig  <- import.bw(con = f, which = nucRef)
    hits <- findOverlaps(nucRef, sig)
    q    <- queryHits(hits)
    s    <- subjectHits(hits)
    st   <- pmax(start(sig)[s], start(nucRef)[q])
    en   <- pmin(end(sig)[s],   end(nucRef)[q])
    res  <- list(pos    = (en + st) / 2,
                 score  = sig$score[s],
                 offset = c(0, cumsum(as.numeric(tabulate(q, nbins = n_nuc)))))
    rm(sig, hits, q, s, st, en)
    gc(verbose = FALSE)
    res
})


## summit position and score of one nucleosome in one sample (NA if none)
get_summit <- function(pos, score){
    # in case no signal at position available assign NA
    if (!(length(score) > 0)) return(c(NA, NA))

    lm<-loess(score~pos, span=0.6)
    lm.fkt<-predict(lm)

    # get local maxima
    maxima<-pos[which(diff(sign(diff(lm.fkt)))==-2)+1]

    #if more maxima available
    if(length(maxima)>1){
        maxima<-maxima[which.max(lm.fkt[which(diff(sign(diff(lm.fkt)))==-2)+1])]
    }

    # in case no local maxima available assign NA
    if (!(length(maxima) > 0)) return(c(NA, NA))

    # summit position and nuc occupancy at summit
    c(maxima, max(lm.fkt[which(diff(sign(diff(lm.fkt)))==-2)+1]))
}


## all nucleosomes, in parallel blocks; each block returns a [nucs x samples x 2] array
blocks <- split(seq_len(n_nuc), cut(seq_len(n_nuc), breaks = min(n_nuc, workers * 4), labels = FALSE))

summit_block <- function(idx){
    res <- array(NA_real_, dim = c(length(idx), length(files.bw), 2))
    for (k in seq_along(idx)) {
        i <- idx[k]
        for (j in seq_along(files.bw)) {
            sp  <- signal_per_nuc[[j]]
            rec <- seq.int(sp$offset[i] + 1, length.out = sp$offset[i + 1] - sp$offset[i])
            res[k, j, ] <- get_summit(sp$pos[rec], sp$score[rec])
        }
    }
    res
}

block_res <- mclapply(blocks, summit_block, mc.cores = workers, mc.preschedule = FALSE)

# a worker killed from outside (usually the OOM killer) leaves NULL, not a try-error.
# Exit with 137 (SIGKILL) so that the pipeline's errorStrategy retries the task
# with more memory, instead of failing later with an unrelated indexing error.
killed <- vapply(block_res, is.null, logical(1))
if (any(killed)) {
    message("summit detection: ", sum(killed), " of ", length(killed),
            " worker blocks returned no result (worker killed, most likely out of memory)")
    quit(save = "no", status = 137)
}
failed <- vapply(block_res, inherits, logical(1), what = "try-error")
if (any(failed)) stop("summit detection failed: ", as.character(block_res[[which(failed)[1]]]))

# preallocated result matrices (rows = nucleosomes, columns = samples)
max.pos.mx   <- matrix(NA_real_, nrow = n_nuc, ncol = length(files.bw), dimnames = list(NULL, names.bw))
max.score.mx <- max.pos.mx
for (b in seq_along(blocks)) {
    max.pos.mx[blocks[[b]], ]   <- block_res[[b]][, , 1]
    max.score.mx[blocks[[b]], ] <- block_res[[b]][, , 2]
}


#index rownames
rownames(max.score.mx)<-paste(seqnames(nucRef),
                              ranges(nucRef), sep = ":")
rownames(max.pos.mx)<-paste(seqnames(nucRef),
                            ranges(nucRef), sep = ":")

# order colum names
time<-as.numeric(substring(colnames(max.score.mx),2))
ord<-order(time)

max.pos.mx<-max.pos.mx[,ord]
max.score.mx<-max.score.mx[,ord]


#calculate the variance of each peak summit over all timepoints
var_pos<-apply(max.pos.mx,1,function(x) var(x,na.rm = T))
var_pos_sort<-sort(var_pos)


############################# LOESS fit ######
######### normalize ##########

## Rank and normalize variation of nucleosome positions
pos_rank<-rank(var_pos_sort, ties.method = "first")/length(var_pos_sort)
var_pos_norm<-(var_pos_sort/max(var_pos_sort))

# fit the data
#define to take 1000 nucleosomes as window for loess
span.loess.pre<-1000/length(var_pos_norm)
span.loess<-ifelse(span.loess.pre<0.05, span.loess.pre, 0.05)

lm<-loess(var_pos_norm~pos_rank,span=span.loess)
#Curve fitting
loess.line<-predict(lm)

#calculate the first derivative
loess.f1<-diff(loess.line)/diff(pos_rank)
#define a slope cutoff of 3
cutOff<-max(which(loess.f1<3))+1

png("dynSHIFT_selection.png",
        width = 1280, height = 1280, res =300 )
    plot(pos_rank,var_pos_norm,
         ylab="variance score of reference position", 
         xlab = "normalized rank", bty="n", pch=19)
    points(pos_rank[cutOff:length(pos_rank)],
           var_pos_norm[cutOff:length(pos_rank)],
           pch=19, col="aquamarine3")
    lines(pos_rank,loess.line, col="red")
    abline(v=pos_rank[cutOff], col="grey", lty=2)
dev.off()


### select the peaks with highest variance
# nucleosomes without a summit in at least two samples have var NA: sort() drops
# them from var_pos_sort, so drop them from the selection as well (else they end
# up as all-NA rows in positions_varSHIFTS.bed)
idx_var<-var_pos>=var_pos_sort[cutOff] & !is.na(var_pos)

############################################



#get chromosome
chr<-sapply(strsplit(names(var_pos),split = ":"),function(x) x[1])
ranges<-sapply(strsplit(names(var_pos),split = ":"),function(x) x[2])
end<-sapply(strsplit(ranges,split = "-"),function(x) x[2])
start<-sapply(strsplit(ranges,split = "-"),function(x) x[1])

#export bed file and define score by variance
var_peaks_bed<-data.frame(chr=chr,
                            start=start,
                            end=end,
                            name=names(var_pos),
                            varSHIFT=var_pos,
                            strand=".")

write.table(var_peaks_bed[order(var_peaks_bed$varSHIFT,decreasing = T),],
            file="shift_result_table.tsv",
            row.names = FALSE, sep="\t", quote=FALSE,col.names = T)


## export selected regions
shift_peaks_bed<-var_peaks_bed[idx_var,]

## export nucleosome positions
write.table(shift_peaks_bed[order(shift_peaks_bed$varSHIFT,decreasing = T),],
            file="positions_varSHIFTS.bed",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = FALSE)



 
