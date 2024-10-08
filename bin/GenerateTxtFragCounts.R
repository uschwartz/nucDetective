#!/usr/bin/env Rscript

### Package loading ###
library(stringr)

#### Function definition ####
######## DATA LOADING ######
read_files <- function(){
    fastqc_files <- list.files(path = ".", pattern = "_fastqc.zip", 
                        all.files = TRUE, full.names = TRUE, recursive = TRUE)
    fastqc_file_raw<- fastqc_files[-grep("_trimmed_", fastqc_files)]
    fastqc_file_trimmed<- grep("_trimmed_", fastqc_files, value = T)
    alignment_file <- list.files(path = ".", pattern = "_alignment_stats.txt",
                            all.files = TRUE, full.names = TRUE)
    qualimap_file <- list.files(path = ".", pattern = "genome_results.txt", 
                            full.names = TRUE, recursive = TRUE)
    filtered_align <- list.files(path = ".", pattern = "_FiltLog.txt",
                                 all.files = TRUE, full.names = TRUE)
    files <- c(fastqc_file_raw[1],fastqc_file_trimmed[1],
               alignment_file, qualimap_file, filtered_align)
    return(files)
}

read_fastqc <- function(fastqc_file){
    fastqc_vec <- c()
    fast_list <- unzip(fastqc_file, list = TRUE)
    for (i in 1:length(stringr::str_detect(fast_list$Name, "/fastqc_data.txt"))){
      if (stringr::str_detect(fast_list$Name, "/fastqc_data.txt")[i] == TRUE) {
        fast_name <- c(fast_list[i,1]) 
      }
    }
    fast <- unzip(fastqc_file, files = fast_name)
    df <- read.table(fast[1], header = FALSE, nrows = 7, fill = TRUE)
    fastqc_vec <- strtoi(df[5,3])
    fastqc_matrix <- matrix(fastqc_vec, dimnames = list("", "Sequenced"))
    return(fastqc_matrix)
}

read_fastqc_trim <- function(fastqc_file, read_frame){
    fastqc_vec <- c()
    fast_list <- unzip(fastqc_file, list = TRUE)
    for (i in 1:length(stringr::str_detect(fast_list$Name, "/fastqc_data.txt"))){
        if (stringr::str_detect(fast_list$Name, "/fastqc_data.txt")[i] == TRUE) {
            fast_name <- c(fast_list[i,1]) 
        }
    }
    fast <- unzip(fastqc_file, files = fast_name)
    df <- read.table(fast[1], header = FALSE, nrows = 7, fill = TRUE)
    fastqc_vec <- strtoi(df[5,3])
    fastqc_matrix <- matrix(fastqc_vec, dimnames = list("", "Trimmed"))
    return_frame <- cbind(read_frame, fastqc_matrix)
    return(return_frame)
}

read_alignedqc <- function(alignment_file, read_frame){
    aligned_reads <- c()
    df <- read.table(alignment_file, header = FALSE, skip = 4, fill = TRUE)
    aligned_reads <- sum(strtoi(df[4,1]), strtoi(df[5,1]))
    return_frame <- data.frame(aligned_reads)
    colnames(return_frame) <- c("Aligned")
    sample <- stringr::str_remove(alignment_file, pattern ="_alignment_stats.txt")
    sample <- stringr::str_remove(sample, pattern = "./")
    rownames(return_frame) <- sample
    return_frame <- cbind(read_frame, return_frame)
    return(return_frame)
}

read_qualimap<- function(qualimap_file, read_frame){
    qual_mat <- c()
    df <- read.table(qualimap_file, header = FALSE, skip = 14, nrows = 7, fill = TRUE)
    qual_vec <- strtoi(df[3,5] %>% stringr::str_remove_all(","))/2
    qual_mat <- rbind(qual_mat, qual_vec)
    colnames(qual_mat) <- c("MAPQC-filtered")
    return_frame <- cbind(read_frame, qual_mat)
    return(return_frame)
}

read_filtered_alignments <- function(f_align_files, read_frame){
    remaining_reads <- matrix(data = NA, nrow = 1, ncol = length(f_align_files), dimnames = list(" ", c( "Nucleosomes")))
    for (j in 1:length(f_align_files)){
      df <- read.delim(f_align_files[j], header = TRUE, skip = 1)
      remaining_reads[1, j] <- df$Reads.Remaining/2
    }
    return_frame <- cbind(read_frame, remaining_reads)
    return_frame <- as.data.frame(return_frame)
    return(return_frame)
}

plot_frag_abs <- function(read_frame){
  plot_name <- paste(row.names(read_frame), c("fragments.pdf"), sep = "_")
  color <- "coral3"
  pdf(plot_name, width =8)
    barplot(as.matrix(read_frame/10**6), beside = FALSE, names.arg = colnames(read_frame), sub = paste("Sample:", row.names(read_frame), sep = " "), 
            col = color, border = NA, ylab = "Number of fragments per million", col.axis = "grey30", col.lab = "grey30",
            col.sub = "grey30", las = 1, ylim = c(0, signif(max(read_frame/10**6), digits = -1)))
  dev.off()
}

### Output Generation ###
files <- read_files()
reads <- read_fastqc(files[1])
reads <- read_fastqc_trim(files[2],reads)
reads <- read_alignedqc(files[3], reads)
reads <- read_qualimap(files[4], reads)
reads <- read_filtered_alignments(files[5], reads)
plot_frag_abs(reads)

txt_name <- rownames(reads)
txt <- write.table(x = reads, 
                   file = paste(txt_name, c("fragment_statistic.txt"), sep = "_") , 
                   append = FALSE, quote = FALSE, sep = "\t", dec = ".",
                   row.names = TRUE, col.names = TRUE)
