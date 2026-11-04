process reference_map{
  container 'uschwartz/r_core:v4.2'
  // resources: see process withName:'reference_map|reference_map_all' in nextflow.config
  publishDir "${params.outDir}/RUN/03_NUCS/referenceMap/", mode: 'copy'

  input:
  file('*')
  val(sampleNames)

  output:
  file ("NucPosRef_top20.bed")

  script:
  """
  NucReferenceMaps.R $params.nucref_min_support $sampleNames
  """
}

process reference_map_all{
  container 'uschwartz/r_nucmacc:v3.1'
  // resources: see process withName:'reference_map|reference_map_all' in nextflow.config
  publishDir "${params.outDir}/RUN/03_NUCS/referenceMap/", mode: 'copy'

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
