process score_peaks{
  label 'mid'
  publishDir "${params.outDir}/RUN/04_NUC_SCORES/", mode: 'copy'

  input:
  file(bw)
  file(nucPos)

  output:
  file("top20_Nucs_Scores.txt")
  file("std.npz")

  script:
  """
  multiBigwigSummary BED-file -b $bw -o std.npz \
  --BED $nucPos --smartLabels \
  --outRawCounts top20_Nucs_Scores.txt -p $task.cpus
  """
}

process pca{
  label 'mid'
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
