process qualimap {

  container 'uschwartz/qualimap'
  label 'mid'
  // resources: see process withName:'qualimap' in nextflow.config
  publishDir "${params.outDir}/QC/05_QUALIMAP", mode: 'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  file "${sampleID}"
  tuple val(sampleID), file("${sampleID}/*.txt")

  script:
  """
  qualimap bamqc --java-mem-size=${task.memory.toGiga() - 2}G -bam $bam -c -outdir ${sampleID}
  """
}
