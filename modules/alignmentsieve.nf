process sieve{
  label 'big'
  memory { params.genomeSize > 200000000 ? '15.GB' : '7.GB'}
  publishDir "${params.outDir}/QC/07_ALIGNMENT_FILTERING/", mode: 'copy', pattern: "*_FiltLog.txt"
  publishDir "${params.outDir}/RUN/00_ALIGNMENT/monoNuc", mode: 'copy', pattern: "*_monoNuc.bam", enabled:params.publishBamFlt

  input:
  tuple val(sampleID), val(bam)

  output:
  tuple val(sampleID), file("*_FiltLog.txt")
  tuple val(sampleID), file("*_monoNuc.bam")

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
  """
}
