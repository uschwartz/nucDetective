#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)
args <- commandArgs(TRUE)

values_files <- list.files()

values2plot<- read.csv(values_files[1])
for (i in 2:length(values_files)){
  values2plot <- rbind(values2plot, read.csv(values_files[i]))
}

write.csv(values2plot, paste0("nrl_values_all.csv"), quote = F, row.names = F)

ylim1<- boxplot.stats(values2plot$nrls)$stats[c(1,5)]
meds <- ddply(values2plot, .(condition), summarise, med = round(median(nrls, na.rm = T),digits=2))

pdf("Comparison_of_NRLs_var_cond.pdf")
ggplot(values2plot, aes(x = condition, y = nrls)) +
  geom_boxplot()+
  coord_cartesian(ylim = ylim1*1.05)+
  labs(title = paste0("Comparison of NRLs under different conditions"))+ theme_bw()+ ylab("NRLs")+
  geom_text(data = meds, aes(x = condition, y = med, label = med), 
            size = 3, vjust = -1.4)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
