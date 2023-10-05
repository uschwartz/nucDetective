#!/usr/bin/env Rscript
library(rtracklayer)

# load reference nucleosomes top20
nucRef<-import.bed("NucPosRef_top20.bed")


files.bw<-grep("bw",list.files(), value = T)
#get name
names.bw<-gsub("_qNorm.bw","",files.bw)


for (i in 1:length(nucRef)){
    #i=1
    #print(i)
    
    #select peak
    queryRef<-nucRef[i]
    
    #get nuc occupancy values from bigwigs
    
    for(j in 1:length(names.bw)){
        bw.vals<-import.bw(con=files.bw[j],
                           which=queryRef)
        #############
        #get midpoint of data
        pos<-(end(bw.vals)+start(bw.vals))/2
        lm<-loess(bw.vals$score~pos, span=0.6)
        lm.fkt<-predict(lm)
        
        # get local maxima
        maxima<-pos[which(diff(sign(diff(lm.fkt)))==-2)+1]
        
        #if more maxima available
        if(length(maxima)>1){
            maxima<-maxima[which.max(lm.fkt[which(diff(sign(diff(lm.fkt)))==-2)+1])]
        }
        
        # in case no local maxima available assign NA
        if(!(length(maxima)>0)){
            max.pos<-NA
            max.score<-NA
            
            ##store results in vector
            if(j==1) max.pos.v<-max.pos else max.pos.v<-c(max.pos.v,max.pos)
            if(j==1) max.score.v<-max.score else max.score.v<-c(max.score.v,max.score)
        } else {
            # nuc occupancy at summit
            max.score<-max(lm.fkt[which(diff(sign(diff(lm.fkt)))==-2)+1])
            
            #select summit position
            max.pos<-maxima
            
            ##store results in vector
            if(j==1) max.pos.v<-max.pos else max.pos.v<-c(max.pos.v,max.pos)
            if(j==1) max.score.v<-max.score else max.score.v<-c(max.score.v,max.score)
        }
        

    }
    # assign timepoints to vectors
    names(max.score.v)<-names.bw
    names(max.pos.v)<-names.bw
    
    #store to matrix
    if(i==1) {
        max.pos.mx<-max.pos.v
        max.score.mx<-max.score.v 
    } else {
        max.pos.mx<-rbind(max.pos.mx, max.pos.v)
        max.score.mx<-rbind(max.score.mx, max.score.v)
    }
    
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
#define a slope cutoff of 1
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
idx_var<-var_pos>=var_pos_sort[cutOff]

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



 
