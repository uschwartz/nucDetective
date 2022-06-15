process nucs2bed{
  container 'uschwartz/r_nucmacc:v3.1'
  publishDir "${params.outDir}/RUN/02_NUCS/", mode: 'copy'

  input:
  tuple val(sampleID), file(nucs)

  output:
  file("**_top20.bed")
  file("bed_allNucs/*.bed")

  script:
  """
  mkdir bed_allNucs
  mkdir bed_top20
  nucs2bedFormat.R $sampleID
  """
}
