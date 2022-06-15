#!/usr/bin/env Rscript

args   <- commandArgs(TRUE)
samp.name<-args[1]

## show DANPOS output files
nuc.xls<-grep(".xls",list.files(),value = T)


for(i in nuc.xls){
    #read DANPOS xls file
    data.xls<-read.delim(i)
    #convert to bed
    data.bed<-data.frame(data.xls[,1:3],
                         summitValue_Position=paste(round(data.xls$smt_value,
                                                          digits = 2),
                                                    data.xls$smt_pos, sep="_"),
                         fuzziness=data.xls$fuzziness_score,
                         strand=rep(".",length(data.xls$chr)))
    
    #write bed
    write.table(data.bed, 
                file =paste0("bed_allNucs/", samp.name,".bed"), 
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote=FALSE )
    
    #get best 20% positioned nucleosomes 
    numNucs<-round(nrow(data.xls)*0.20)
    
    #write bed
    write.table(data.bed[1:numNucs,], 
                file =paste0("bed_top20/", samp.name,"_top20.bed"), 
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote=FALSE )
    
} 
