#!/usr/bin/env Rscript
# Relative occurrence of dynamic nucleosome features (regularity changes, occupancy,
# shifts, fuzziness) along scaled gene bodies.
#
# usage: plot_feature_distribution.R <values_featureDistribution.txt> [upstream=750] [downstream=750]
# input: plotProfile --outFileNameData table of a computeMatrix scale-regions matrix

suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
})

args <- commandArgs(TRUE)
infile <- args[1]
upstream <- if (length(args) >= 2) as.numeric(args[2]) else 750
downstream <- if (length(args) >= 3) as.numeric(args[3]) else 750

lines <- strsplit(readLines(infile), "\t")

# row 1: bin labels (tick labels such as TSS / TES), row 2: bin numbers,
# following rows: <feature> <group> <values per bin>.
# With --perGroup, deepTools repeats the label block once per sample, so only the
# first block (as long as one data row) is used.
rows <- lines[-(1:2)]
features <- vapply(rows, `[`, "", 1)
mat <- do.call(rbind, lapply(rows, function(r) as.numeric(r[-(1:2)])))
rownames(mat) <- features
n.bins <- ncol(mat)
bin.labels <- c(lines[[1]][-(1:2)], rep("", n.bins))[1:n.bins]

tss.bin <- which(bin.labels == "TSS")
tes.bin <- which(bin.labels == "TES")

# relative occurrence: each feature sums up to 1 over the whole profile
total <- rowSums(mat)
empty <- !(total > 0)
if (any(empty)) {
    message("No signal over genes for: ", paste(features[empty], collapse = ", "), " (not plotted)")
}
rel <- mat[!empty, , drop = FALSE] / total[!empty]

df <- melt(rel, varnames = c("feature", "position"), value.name = "relative_occurrence")
df$feature <- as.character(df$feature)
write.table(df, "relative_featureDistribution.txt", sep = "\t", quote = FALSE, row.names = FALSE)

feature.colors <- c(fuzziness = "brown2", occupancy = "orange",
                    regularity_changes = "#7F81B0", shifts = "aquamarine3")

p <- ggplot(df, aes(x = position, y = relative_occurrence, color = feature)) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.1, se = FALSE) +
    theme_classic() +
    scale_color_manual(values = feature.colors) +
    ylab("relative occurrence rate") +
    xlab("") +
    scale_x_continuous(breaks = c(1, tss.bin, tes.bin, n.bins),
                       labels = c(-upstream, "TSS", "TES", paste0("+", downstream))) +
    geom_vline(xintercept = c(tss.bin, tes.bin), color = "gray", linetype = "dashed")

pdf("featureDistribution.pdf", width = 5, height = 3)
print(p)
invisible(dev.off())
