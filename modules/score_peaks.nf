process score_peaks{
  label 'big'
  // resources: see process withName:'score_peaks' in nextflow.config
  // chrSizes is read (see script: out-of-bounds clip) but passed as a param path,
  // so bind-mount its directory - same pattern as regularity.nf.
  if(params.container_engine == 'docker'){
    containerOptions "-v \$(dirname ${params.chrSizes}):\$(dirname ${params.chrSizes})"
  }
  if(params.container_engine == 'singularity'){
    containerOptions "-B \$(dirname ${params.chrSizes}):\$(dirname ${params.chrSizes})"
  }
  publishDir "${params.outDir}/RUN/04_NUC_SCORES/", mode: 'copy'

  input:
  file(bw)
  file(nucPos)

  output:
  file("top20_Nucs_Scores.txt")
  file("std.npz")

  script:
  """
  # multiBigwigSummary BED-file's -p only splits work by chromosome, so on a
  # large reference map one process grinds through all of chr1 while the rest
  # sit idle (observed: ~1 core busy for a day on a ~2.4 M-nucleosome map).
  # Split the BED into \$task.cpus equal parts, run one single-threaded
  # multiBigwigSummary per part in parallel, then stitch the parts back
  # together. Output is identical: multiBigwigSummary keeps BED input order and
  # the parts are concatenated in split order; the per-part .npz matrices are
  # row-stacked (deeptools npz = matrix + labels), which plotPCA / plotCorrelation
  # read unchanged.
  N=${task.cpus}
  mkdir -p _parts

  # A few reference nucleosomes have their +-50 bp extension (NucReferenceMaps.R)
  # running past a scaffold end; multiBigwigSummary aborts on those with
  # "Invalid interval bounds!". Clip every region's end to the chromosome length
  # (17 / 2.4 M on GRCh38, all on scaffolds far from the cell-line signal).
  awk 'BEGIN{OFS="\\t"}
       NR==FNR { L[\$1]=\$2; next }
       (\$1 in L) { if (\$3 > L[\$1]) \$3 = L[\$1]; if (\$2 < \$3) print }' \
      $params.chrSizes $nucPos > nucPos.clipped.bed

  split -d -a 4 --additional-suffix=.bed -n l/\$N nucPos.clipped.bed _parts/p

  pids=""
  for part in _parts/p*.bed; do
    base="\${part%.bed}"
    multiBigwigSummary BED-file -b $bw -o "\${base}.npz" \
      --BED "\$part" --smartLabels --outRawCounts "\${base}.tab" -p 1 &
    pids="\$pids \$!"
  done
  for p in \$pids; do wait "\$p"; done

  tabs=(_parts/p*.tab)
  head -n 1 "\${tabs[0]}" > top20_Nucs_Scores.txt
  for f in "\${tabs[@]}"; do tail -n +2 "\$f" >> top20_Nucs_Scores.txt; done

  merge_npz.py '_parts/p*.npz' std.npz

  rm -rf _parts
  """
}

process pca{
  label 'mid'
  // resources: see process withName:'pca|correlation' in nextflow.config
  publishDir "${params.outDir}/RUN/05_NUC_Explorative/", mode: 'copy'

  input:
  file(scores)


  output:
  file ("PCA.pdf")

  script:
  """
  plotPCA --transpose --corData $scores \
  -o PCA.pdf
  """
}

process correlation{
  label 'mid'
  // resources: see process withName:'pca|correlation' in nextflow.config
  publishDir "${params.outDir}/RUN/05_NUC_Explorative/", mode: 'copy'

  input:
  file(scores)


  output:
  file ("corHeat.pdf")

  script:
  """
  plotCorrelation --corMethod spearman --corData $scores \
  --whatToPlot heatmap  -o corHeat.pdf
  """
}
