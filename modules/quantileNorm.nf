process quantNorm{
  container 'uschwartz/danpos'
  memory { params.genomeSize > 200000000 ? '47.GB' : '15.GB' }
  publishDir "${params.outDir}/RUN/01_NORM_PROFILE", mode: 'copy',
  pattern: "*_qNorm.bw"

  input:
  tuple val(sampleID), file(wig), val(type), val(refID), file(ref), val(ref_type)


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

process wig_to_bw{
  container 'uschwartz/danpos'
  publishDir "${params.outDir}/RUN/01_BIGWIG_PROFILES", mode: 'copy', pattern: "*.bw"

  input:
  tuple val(sampleID), file(wig)
  output:
  tuple val(sampleID), file("*.bw")


  script:
  """
  wigToBigWig $wig -clip $params.chrSizes $sampleID".bw"
  """
}
