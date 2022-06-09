process InsertSize_Histogram{
  container 'uschwartz/r_nucmacc:v3.1'

  publishDir "${params.outDir}/QC/04_FRAGMENT_SIZES", mode: 'copy'

  input:
  file('*')

  output:
  file("*.pdf")

  script:
  """
  InsertSize_Histogram.R $params.minLen $params.maxLen
  """
}
