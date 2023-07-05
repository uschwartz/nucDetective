process qualimap {

  container 'uschwartz/qualimap'
  label 'mid'
  memory '7.GB'
  publishDir "${params.outDir}/QC/05_QUALIMAP", mode: 'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  file "${sampleID}"
  tuple val(sampleID), file("${sampleID}/*.txt")

  script:
  """
  qualimap bamqc --java-mem-size=7G -bam $bam -c -outdir ${sampleID}
  """
}
