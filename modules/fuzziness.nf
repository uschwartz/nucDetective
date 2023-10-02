process fuzziness{
  container 'leoschmutterer/fuzziness:v1.0'
  memory { params.genomeSize > 200000000 ? '60.GB' : '30.GB'}
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
