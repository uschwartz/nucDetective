process quantNorm{
  container 'uschwartz/danpos'
  memory { params.genomeSize > 200000000 ? 50.GB : 16.GB }
  publishDir "${params.outDir}/RUN/01_NORM_PROFILE", mode: 'copy',
  pattern: "*_qNorm.bw"

  input:
  tuple val(sampleID), file(wig), val(type)
  tuple val(refID), file(ref), val(ref_type)


  output:
  tuple val(sampleID), file("*_qNorm.bw")
  tuple val(sampleID), file("wiq_result/*.wig")

  script:
  """
  danpos.py wiq $params.chrSizes $wig --reference $ref
  wigToBigWig wiq_result/*.wig -clip $params.chrSizes $sampleID"_qNorm.bw"
  """
}


process ref_bw{
  container 'uschwartz/danpos'
  publishDir "${params.outDir}/RUN/01_NORM_PROFILE", mode: 'copy',
  pattern: "*_qNorm.bw"

  input:
  tuple val(sampleID), file(ref), val(type)
  output:
  tuple val(sampleID), file("*_qNorm.bw")
  tuple val(sampleID), file("*_qNorm.wig")


  script:
  """
  wigToBigWig $ref -clip $params.chrSizes $sampleID"_qNorm.bw"
  mv $ref $sampleID"_qNorm.wig"
  """
}
