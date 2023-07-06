
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


# Main --------------------------------------------------------------------

# Load Reference
refPos <- read_bed("NucPosRef_top20.bed")

# Load Nucleosome and Fuzziness data
inFiles <- grep("bed",list.files(), value = T)
fileNames <- c("T10", "T15", "T20", "T25", "T30", "T35", "T40", "T5")


inDat <- inFiles %>%
  set_names(fileNames) %>%
  map_df(~read_bed(.x) %>% 
           as_tibble(), .id = "timepoints") %>% 
  separate(name, into = c("occ","dyad"), sep = "_") %>% 
  mutate(occ = as.numeric(occ), fuz = score)

# Select timepoint and fuzziness score from nucleosome data
fuz2 <- inDat %>% 
  as_granges() %>%
  plyranges::select(fuz, timepoints)

# Get overlap of reference with timepoint nucleosome data
hits <- findOverlaps(refPos,fuz2)
gr.over <- pintersect(refPos[queryHits(hits)],fuz2[subjectHits(hits)])
w <- width(gr.over)

# Join reference with nucleosome fuzziness data and select for each reference 
# position the nucleosome with the largest overlap
refFuz <- refPos %>% 
  mutate(nucID = paste0("nuc",1:length(refPos))) %>% 
  join_overlap_inner(fuz2) %>%
  mutate(overlap = w) %>%
  as_tibble() %>%
  group_by(nucID, timepoints) %>% 
  filter(overlap == max(overlap)) %>% 
  # dplyr::select(nucID,fuz, timepoints) %>% 
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
slope_at_cutoff = 1
fitx <- smooth.spline(x = ranked_dat$var_rank, y = ranked_dat$norm_var)
f <-predict(fitx)
f1<-predict(fitx, deriv=1)
# cutoff is the position where the slope is > slope_at_cutoff first time after the 1000 lowest variances. 
# This is necessary as sometimes the first ones have a quite high predicted slope. 
cutOff<-(which(f1$y>slope_at_cutoff)[min(which(which(f1$y>slope_at_cutoff) > 1000))])/nrow(ranked_dat)


# Extract nucIDs of nucleosomes with high variance in fuzziness
highFuzNucs <- ranked_dat %>% filter(var_rank>cutOff) %>% pull(nucID)


# Plot 
ranked_dat %>% ggplot(aes(x = var_rank, y = norm_var)) +
  geom_point() +
  geom_point(data = ranked_dat %>% filter(var_rank>cutOff), aes(x = var_rank, y = norm_var), color = "#EF6461") +
#  geom_smooth(aes(x = var_rank, y = norm_var), method = "loess", span = 0.01, color = "black", size = 0.5) +
  geom_vline(xintercept = cutOff, linetype = "dashed") +
  theme_classic() +
  labs(x = "normalized rank", 
       y = "norm. variance of fuzziness scores") +
  labs_pubr()

ggsave(paste0("var_cutoff_DANPOS.pdf"), width = 4, height = 4)


# Join Data into one table and clean
refFuzTableAllTP <- refFuz %>% 
  left_join(ranked_dat) %>%
  mutate(fuzzinessVariance = var,
         fuzzinessScore = score, 
         highFuzVariance = if_else((nucID %in% highFuzNucs), TRUE, FALSE)) %>% 
  ungroup() %>%
  dplyr::select(-c(overlap, fuz, var_rank, norm_var, var))

# Save Data table with all relevant datacolumns
refFuzTableAllTP %>% 
  dplyr::select(-c(width, strand, name, score)) %>% 
  write_delim(file = paste0("referenceMapped_Fuzziness_table.tsv", delim = "\t"))

# Save Bed file with high variance in fuzziness nucleosome positions
refFuzTableAllTP %>% 
  dplyr::select(-c(name, score, timepoints, fuzzinessScore)) %>% 
  distinct() %>% 
  mutate(name = nucID,
         score = fuzzinessVariance) %>% 
  filter(highFuzVariance) %>% 
  as_granges() %>% 
  write_bed(file = paste0("/HighFuzNuc.bed"))
