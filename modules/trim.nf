process trimgalore{
  label 'big'
  container 'uschwartz/core_docker:v1.0'
  publishDir "${params.outDir}/QC/02_TRIMMING/${sampleID}", mode: 'copy', pattern: "*trimming_report.txt"
  publishDir "${params.outDir}/QC/03_TRIMMED_FASTQC/${sampleID}", mode: 'copy', pattern: "*.html"


  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
  tuple val(sampleID), file("*_1.fq.gz"), file("*_2.fq.gz")
  file "*trimming_report.txt"
  file "*_fastqc.zip"
  file "*_fastqc.html"
  tuple val(sampleID), file("*_fastqc.zip")


  script:
  """
  trim_galore --gzip \
   --paired \
   --cores $task.cpus \
   --fastqc \
   -q 10 \
   --stringency 2 \
   $read1 $read2

   rename 's/_val_/_trimmed_/' *_val_*
  """

}
