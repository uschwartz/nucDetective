#!/usr/bin/env bash
# Distribution of dynamic nucleosome features along genes (scale-regions).
#
# usage: feature_distribution_matrix.sh <annotation.gtf> <chrom_sizes> <cpus> label1:file1.bed [label2:file2.bed ...]
#
# Each feature BED is cleaned (rows with non-numeric coordinates, e.g. "NA", are dropped),
# merged and converted to a 0/1 coverage bigWig. The bigWigs are then summed over genes
# with deepTools computeMatrix scale-regions (-b 750 -a 750) and plotProfile (--averageType sum).
# Output: values_featureDistribution.txt (+ deepTools matrix and default png)
set -euo pipefail

gtf=$1
chrSizes=$2
cpus=$3
shift 3

bws=()
labels=()
for entry in "$@"; do
    label=${entry%%:*}
    bed=${entry#*:}

    awk 'BEGIN{OFS="\t"} $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $1,$2,$3}' "$bed" \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        | bedtools genomecov -bga -i - -g "$chrSizes" \
        | sort -k1,1 -k2,2n > "${label}.bedgraph"

    wigToBigWig -clip "${label}.bedgraph" "$chrSizes" "${label}.bw"

    bws+=("${label}.bw")
    labels+=("$label")
done

computeMatrix scale-regions -S "${bws[@]}" \
    -R "$gtf" \
    -o computeMatrix_featureDistribution.txt.gz \
    --averageTypeBins sum \
    -b 750 -a 750 \
    --samplesLabel "${labels[@]}" \
    -p "$cpus"

plotProfile -m computeMatrix_featureDistribution.txt.gz \
    -out DefaultProfile_featureDistribution.png \
    --averageType sum \
    --outFileNameData values_featureDistribution.txt \
    --perGroup
