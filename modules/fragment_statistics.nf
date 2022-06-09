process statistics_read{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/QC/08_FRAGMENT_STATISTICS/${sampleID}", mode: 'copy'

  input:
  tuple val(sampleID), file(sieve_mono), file(fastqc), file(trim),
   file(alignment), file(qualimap)

  output:
  file("*.txt")
  file("*.pdf")

  script:
  """
  GenerateTxtFragCounts.R
  """
}

process statistics_plot{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/QC/08_FRAGMENT_STATISTICS", mode: 'copy'

  input:
  file(statistics_read)

  output:
  file("*.txt")
  file("*.pdf")

  script:
  """
  Plot_comparison.R
  """
}
