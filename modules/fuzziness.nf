process fuzziness{
  container 'uschwartz/r_core:v1.0'
  publishDir "${params.outDir}/RUN/07_FUZZY_NUCS/", mode: 'copy'

  input:
  file(nucPos)
  file(bed)

  output:
  file ("*")

  script:
  """
  fuzziness.R
  """
}
