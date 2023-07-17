process alignment{

  cpus = (Runtime.runtime.availableProcessors() - (Math.round((Runtime.runtime.availableProcessors()*20)/100)))
  memory = { 2.GB * task.cpus }
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
