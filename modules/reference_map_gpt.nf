process reference_map_gpt{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/RUN/03_NUCS/referenceMap/gpt", mode: 'copy'

  input:
  file('*')
  val(sampleNames)

  output:
  file ("NucPosRef_top20.bed")

  script:
  """
  referencemap_gpt.R $sampleNames
  """
}

process reference_map_all_gpt{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/RUN/03_NUCS/referenceMap/", mode: 'copy'

  input:
  file('*')
  val(sampleNames)

  output:
  file ("NucPosRef_allNucs.bed")

  script:
  """
  referencemap_gpt_all.R $sampleNames
  """
}
