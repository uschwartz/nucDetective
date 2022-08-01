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

# save objects
save(max.score.mx, file="max.score.mx.rda")
save(max.pos.mx, file="max.pos.mx.rda")



