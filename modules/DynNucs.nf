process occupancy{
  container 'uschwartz/r_core:v1.0'
  publishDir "${params.outDir}/RUN/06_DYNAMIC_NUCS/occupancy/", mode: 'copy'

  input:
  file(scores)

  output:
  file ("*")

  script:
  """
  getDynOccNucs.R
  """
}

process shift{
  container 'uschwartz/r_core:v1.0'
  publishDir "${params.outDir}/RUN/06_DYNAMIC_NUCS/shift/", mode: 'copy'

  input:
  file(nucPos)
  file(bw)

  output:
  file ("*")

  script:
  """
  getNucShift.R
  """
}
