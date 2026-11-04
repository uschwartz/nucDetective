process fuzziness{
  container 'leoschmutterer/fuzziness:v1.0'
  // resources: see process withName:'fuzziness' in nextflow.config
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
