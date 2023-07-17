process fuzziness{
  container 'uschwartz/r_core:v4.2'
  publishDir "${params.outDir}/RUN/06_DYNAMIC_NUCS/fuzziness", mode: 'copy'

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
