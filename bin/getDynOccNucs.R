#!/usr/bin/env Rscript

library(DESeq2)
# load score table for reference nucleosomes
peakScores<-read.table("top20_Nucs_Scores.txt",
                       comment.char = "", header = T)
#split position and values
pos<-(peakScores[,1:3])
colnames(pos)<-c("chr","start","end")
pos$strand<-"."


val<-(peakScores[,4:ncol(peakScores)])
colnames(val)<-gsub("_qNorm","",colnames(val))

## Analysis
### get most variable peaks

#rlog transformation
mx<-round(as.matrix(val)*100,digits = 0)
rld<-rlog(mx)

#calculate the variance of each peak over all timepoints
var_counts<-apply(rld,1, var)
names(var_counts)<-paste0("nuc",1:length(var_counts))
var_counts_sort<-sort(var_counts)

# define most dynamic nucleosomes based on slope = 1 cut-off 


##############  LOESS FIT ###################
counts_rank<-rank(var_counts_sort, ties.method = "first")/length(var_counts_sort)
var_counts_norm<-(var_counts_sort/max(var_counts_sort))

### get LOESS smoothing function 
#define to take 1000 nucleosomes as window for loess
span.loess.pre<-1000/length(var_counts_norm)
span.loess<-ifelse(span.loess.pre<0.05, span.loess.pre, 0.05)

lm<-loess(var_counts_norm~counts_rank,span=span.loess)
#Curve fitting
loess.line<-predict(lm)

#get first derivative to deduce the slope of loess function
loess.f1<-diff(loess.line)/diff(counts_rank)
#set the slope cutOff to 1
cutOff<-max(which(loess.f1<3))+1


png("dynOCC_selection.png",
    width = 1280, height = 1280, res =300 )
plot(counts_rank,var_counts_norm,
     ylab="variance score of nucleosome occupancy", 
     xlab = "normalized rank", bty="n", pch=19)
points(counts_rank[cutOff:length(counts_rank)],
       var_counts_norm[cutOff:length(counts_rank)],
       pch=19, col="orange")
lines(counts_rank,loess.line, col="red")
abline(v=counts_rank[cutOff], col="grey", lty=2)
dev.off()


### select the peaks with highest variance
idx_var<-var_counts>=var_counts_sort[cutOff]

### generate data table
df<-data.frame(nucID=names(var_counts), varOCC=var_counts)

df$category<-"normal"
df$category[idx_var]<-"dynamicOccupancy"


##concatenate tables
nucStats<-cbind(pos,df)

nucStats.ord<-nucStats[order(nucStats$varOCC,decreasing = T),]

write.table(nucStats.ord, file="occupancy_result_table.tsv",
            row.names = FALSE, sep="\t", quote=FALSE,col.names = T)


## export nucleosome positions

final.pos<-nucStats.ord[nucStats.ord$category =="dynamicOccupancy",
c("chr", "start","end", "nucID","varOCC", "strand")]

write.table(final.pos[order(final.pos$varOCC,decreasing = T),],
            file="positions_varOCC.bed",
            row.names = FALSE, sep="\t", quote=FALSE, col.names = FALSE)
