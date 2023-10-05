#!/usr/bin/env Rscript
# Fuzziness Metric - DANPOS --------------------------------------------------

# Script by Simon Holzinger mod. by Leo Schmutterer to fit the pipeline
#
# Requires nucleosome position reference and individual nucleosome position files with fuzziness scores
#
# Uses ´slope_at_cutoff´ variable as cutoff for high fuzziness
#
# returns cutoff-plot, a table with fuzziness values, fuzziness variance and 
# affiliation to high fuzziness mapped to the respective reference nucleosome
# positions and a bed file of high fuzziness variation nucleosomes.

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(plyranges)
library(stringr)
library(ggpubr)
library(tidyfast)

# Main --------------------------------------------------------------------

# Load Reference
refPos <- read_bed("NucPosRef_top20.bed")

# Load Nucleosome and Fuzziness data
inFiles <- grep("bed",list.files(), value = T)
inFiles<- inFiles[-grep("NucPosRef", inFiles)]
fileNames <- c()
for(i in 1:length(inFiles)) fileNames[i] <- gsub(".bed","", inFiles[i])
fileNames <- fileNames[order(fileNames)]



inDat <- inFiles %>%
  set_names(fileNames) %>%
  map_df(~read_bed(.x) %>% 
           as_tibble(), .id = "timepoints") %>% 
  dt_separate(name, into = c("occ","dyad"), sep = "_") %>% 
  mutate(occ = as.numeric(occ), fuz = score)

# Select timepoint and fuzziness score from nucleosome data
fuz2 <- inDat %>% 
  as_granges() %>%
  plyranges::select(fuz, timepoints)

rm(inDat)
# Get overlap of reference with timepoint nucleosome data
hits <- findOverlaps(refPos,fuz2)
gr.over <- pintersect(refPos[queryHits(hits)],fuz2[subjectHits(hits)])
w <- width(gr.over)

rm(hits, gr.over)
# Join reference with nucleosome fuzziness data and select for each reference 
# position the nucleosome with the largest overlap
refFuz <- refPos %>% 
  mutate(nucID = paste0("nuc",1:length(refPos))) %>% 
  join_overlap_inner(fuz2) %>%
  mutate(overlap = w) %>%
  as_tibble() %>%
  group_by(nucID, timepoints) %>% 
  filter(overlap == max(overlap)) %>% 
  filter(!is.na(fuz)) %>% 
  mutate(score = fuz,
         name = paste(nucID, timepoints, sep = "_"))

# Calculate variance, rank and normalize both 
ranked_dat <- refFuz %>%
  group_by(nucID) %>%
  summarise(var = var(fuz)) %>% 
  ungroup() %>% 
  mutate(var = ifelse(is.na(var),0,var)) %>%  #filter variance == NA 
  arrange(var) %>% 
  mutate(var_rank = 1:dplyr::n(),
         var_rank = var_rank/max(var_rank),
         norm_var = var/max(var))

# fit curve and get slope cutoff
slope_at_cutoff = 3
span.loess.pre<-1000/length(ranked_dat$norm_var)
span.loess<-ifelse(span.loess.pre<0.05, span.loess.pre, 0.05)

lm<-loess(ranked_dat$norm_var~ranked_dat$var_rank,span=span.loess)
#Curve fitting
loess.line<-predict(lm)

#calculate the first derivative
loess.f1<-diff(loess.line)/diff(ranked_dat$var_rank)

cutOff<-max(which(loess.f1<3))+1

cutOffval <- ranked_dat$var_rank[cutOff]

# Extract nucIDs of nucleosomes with high variance in fuzziness
highFuzNucs <- ranked_dat %>% filter(var_rank>=cutOffval) %>% pull(nucID)


# Plot 
png(paste0("dynFUZZ_selection.png"),width = 1280, height = 1280, res =300)
plot(ranked_dat$var_rank, ranked_dat$norm_var,
     ylab = "norm 10xLog of reference position", 
     xlab = "norm rank", bty = "n", pch = 19)
points(ranked_dat$var_rank[cutOff:length(ranked_dat$norm_var)],
       ranked_dat$norm_var[cutOff:length(ranked_dat$norm_var)],
       pch=19, col="#7F81B0")
lines(ranked_dat$var_rank,loess.line, col="red")
abline(v = ranked_dat$var_rank[cutOff], col = "blue", lty = 2)
dev.off()


# Join Data into one table and clean
refFuzTableAllTP <- refFuz %>% 
  left_join(ranked_dat) %>%
  mutate(varFuzz = var,
         fuzzinessScore = score, 
         category = if_else((nucID %in% highFuzNucs), "dynamicFuzziness", "normal"),
         highFuzVariance = if_else((nucID %in% highFuzNucs), T, F)) %>% 
  ungroup() %>%
  dplyr::select(-c(overlap, fuz, var_rank, norm_var, var))

colnames(refFuzTableAllTP)[1] <- "chr"
# Save Data table with all relevant datacolumns
refFuzTableAllTP %>% 
  dplyr::select(-c(width, name, score, timepoints, fuzzinessScore, highFuzVariance))%>%
  distinct() %>% arrange(desc(varFuzz))%>%
  write_delim(file = paste0("fuzziness_result_table.tsv"), delim = "\t")

colnames(refFuzTableAllTP)[1] <- "seqnames"
# Save Bed file with high variance in fuzziness nucleosome positions
refFuzTableAllTP %>% 
  dplyr::select(-c(name, score, timepoints, fuzzinessScore, category)) %>% 
  distinct() %>% 
  mutate(name = nucID,
         score = varFuzz) %>% 
  filter(highFuzVariance) %>% 
  as_granges() %>% 
  write_bed(file = paste0("positions_varFUZZ.bed"))
