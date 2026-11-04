process alignment{

  if (params.container_engine == 'docker') {
    containerOptions "-v \$(dirname ${params.genomeIdx}):\$(dirname ${params.genomeIdx})"
  }
  if(params.container_engine == 'singularity'){
    containerOptions "-B \$(dirname ${params.genomeIdx}):\$(dirname ${params.genomeIdx})"
  }


  // resources (cpus/memory): see process withName:'alignment' in nextflow.config
  publishDir "${params.outDir}/QC/04_ALIGNMENT", mode: 'copy', pattern: "*_alignment_stats.txt"
  publishDir "${params.outDir}/RUN/00_ALIGNMENT", mode: 'copy', pattern: "*_aligned.bam", enabled:params.publishBam



  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
  file "*_alignment_stats.txt"
  tuple val(sampleID), file("*_aligned.bam") 
  tuple val(sampleID), file("*_alignment_stats.txt")
  tuple val(sampleID), file("*_aligned.bam"), file("*.bai")	
  
  script:
  """
  bowtie2 -t \
  --threads $task.cpus \
  --very-sensitive \
  --no-discordant \
  --no-mix \
  --dovetail \
  -x $params.genomeIdx \
  -1 $read1 \
  -2 $read2 \
  2> ${sampleID}_alignment_stats.txt \
  | samtools view -bS -q $params.MAPQC -f 2 -@ $task.cpus - | samtools sort -@ $task.cpus - > ${sampleID}"_aligned.bam"

  samtools index -@ $task.cpus -b ${sampleID}"_aligned.bam" 
  """

}
