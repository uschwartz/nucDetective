process danpos{
  container 'uschwartz/danpos'
  memory { params.genomeSize > 200000000 ? 50.GB : 16.GB }
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
  resolution = ( params.genomeSize > 200000000 ? '10':'1')
  """
  danpos.py dpos $bam -m 1 --extend 70 -c $params.genomeSize \
  -u 0 -z 20 -a $resolution -e 1 \
  --distance 75 --width 10  > $sampleID"_DANPOS_stats.txt"
  wigToBigWig result/pooled/*.wig -clip $chrSizes $sampleID"_monoNucs_profile.bw"
  """
}


//mv  result/pooled/*.wig ./$sampleID"_profile.wig"
