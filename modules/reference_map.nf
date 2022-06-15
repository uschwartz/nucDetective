process reference_map{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/RUN/02_NUCS/referenceMap/", mode: 'copy'

  input:
  file('*')
  val(sampleNames)

  output:
  file ("NucPosRef_top20.bed")

  script:
  """
  NucReferenceMaps.R $sampleNames
  """
}

process reference_map_all{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/RUN/02_NUCS/referenceMap/", mode: 'copy'

  input:
  file('*')
  val(sampleNames)

  output:
  file ("NucPosRef_allNucs.bed")

  script:
  """
  NucReferenceMaps_all.R $sampleNames
  """
}
