#!/usr/bin/env Rscript

#load libraries
library(plyranges)
library(tidyverse)
library(ggpubr)


#load data
dat <- list.files()
NRLs <- bind_rows(read_csv(dat))
NRLs <- NRLs[order(NRLs$sample),]


#plot NRls

plt <- NRLs %>% 
  ggplot() +
  geom_pointrange(aes(ymin = nrl.2.5, ymax = nrl.97.5, x = sample, y = NRL, color = sample), position = position_dodge(width = 0.4)) +
  #scale_color_manual(values = time.pal) +
  labs(title = "Nucleosome Repeat Length",
       x = "", 
       y = "mean NRL [bp]",
       caption = "95%CI indicated") +
  theme_pubr() +
  guides(colour = "none") +
  labs_pubr() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Compare_NRLs.pdf", plot = plt, height = 5, width = 3.5)

write_csv(NRLs, file = paste0("Compare_NRLs.csv"))