process danpos{
  container 'uschwartz/danpos'
  memory '15.GB'
  publishDir "${params.outDir}/RUN/01_NUCLEOSOME_PROFILE", mode: 'copy',
  pattern: "*_monoNucs_profile.bw"
  publishDir "${params.outDir}/RUN/01_NUCLEOSOME_PROFILE/wig", mode: 'copy',
  pattern: "result/pooled/*.wig",saveAs: {filename -> "${sampleID}_profile.wig"},
  enabled:!params.OmitPublishWig

  input:
  tuple val(sampleID), file(bam)
  file(chrSizes)

  output:
  file("*_monoNucs_profile.bw")
  tuple val(sampleID), file("result/pooled/*.xls")
  file("result/pooled/*.wig")

  script:
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 20 -e 1 \
  --distance 75 --width 10  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_monoNucs_profile.bw"
  """
}
