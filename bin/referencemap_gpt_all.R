#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

library(rtracklayer)
library(GenomicRanges)

# load bed files
bed.list <- lapply(args, function(i) {
  name<-paste0(i,".bed")
  print(name)
  import.bed(name)
})

#start with first entry as query
newquery <- bed.list[[1]]
#remove from list
timepoint_peaks <- bed.list[-1]

p_reduce <- GRanges()
# iterate over all timepoints
for(a in seq_along(timepoint_peaks)) {
  print(c("This is round #", a))

  query1 <- newquery
  subject1 <- timepoint_peaks[[a]]

  # create Hits object, to later exclude them from query/ subject
  temp_overlap <- findOverlaps(query1, subject1, minoverlap = 100)

  # create Pairs object to reduce regions pairwaise
  temp_pairs <- findOverlapPairs(query1, subject1, minoverlap = 100)

  # get regions that are only in query/ subject to include them later
  only_A <- query1[-queryHits(temp_overlap)]
  only_B <- subject1[-subjectHits(temp_overlap)]

  # test exact overlap of regions
  exactOvrl <- temp_pairs@first[temp_pairs@first == temp_pairs@second]
  temp_pairs2 <- temp_pairs[!(temp_pairs@first == temp_pairs@second)]

  # reduce pairwise (get region from query.start to subject.end)
  print("Now reducing pairwise, this might take some time")

  s <- reduce(GRanges(seqnames = seqnames(temp_pairs2@first),
                      ranges = IRanges(start = start(ranges(temp_pairs2@first)),
                                       end = end(ranges(temp_pairs2@second)))))

  p_reduce <- c(p_reduce, s)

  # create query for next comparison step with regions: Overlapping, only in query, only in subject
  newquery <- c(p_reduce, exactOvrl ,only_A, only_B)
  # empty GRanges again for next pairwise reduction
  p_reduce <- GRanges()
}

## merge reference peaks if they share 120 bp
nucRef <- newquery
OvrlpCount <- countOverlaps(nucRef,nucRef, minoverlap = 120)
ovrlps <- which(OvrlpCount > 1)

if(length(ovrlps) > 0){
  merged <- reduce(nucRef[ovrlps])
  newquery <- c(nucRef[-ovrlps], merged)
}

# write final result to bed file
base <- as.data.frame(newquery)
# therefore bring to bed file order
base_w <- data.frame(base[,1:3], base$name, base$score, base$strand)
write.table(base_w, file= "NucPosRef_allNucs.bed",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE)
