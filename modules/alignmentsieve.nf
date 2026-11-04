process sieve{
container 'uschwartz/deeptools_samtools:v1.0'

if(params.blacklist && params.container_engine == 'docker'){
    containerOptions "-v \$(dirname ${params.blacklist}):\$(dirname ${params.blacklist})"
}
if(params.blacklist && params.container_engine == 'singularity'){
    containerOptions "-B \$(dirname ${params.blacklist}):\$(dirname ${params.blacklist})"
}

  label 'big'
  // resources: see process withName:'sieve' in nextflow.config
  publishDir "${params.outDir}/QC/07_ALIGNMENT_FILTERING/", mode: 'copy', pattern: "*_FiltLog.txt"
  publishDir "${params.outDir}/RUN/00_ALIGNMENT/monoNuc", mode: 'copy', pattern: "*_monoNuc.bam", enabled:params.publishBamFlt

  input:
  tuple val(sampleID), file(bam), file(idx)

  output:
  tuple val(sampleID), file("*_FiltLog.txt")
  tuple val(sampleID), file("*_monoNuc.bam")
  tuple val(sampleID), file("*_monoNuc.bam"), file("*.bai")

  script:
  blacklistOpt = ( params.blacklist ? "--blackListFileName $params.blacklist":'')
  """
  alignmentSieve -b $bam \
  -o ${sampleID}"_monoNuc.bam" \
  -p $task.cpus \
  --filterMetrics  ${sampleID}"_FiltLog.txt" \
  --minFragmentLength $params.minLen \
  --maxFragmentLength $params.maxLen \
  $blacklistOpt
  
  samtools index -@ $task.cpus -b ${sampleID}"_monoNuc.bam"
  """
}
