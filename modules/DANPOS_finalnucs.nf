process danpos_nuc{
  container 'uschwartz/danpos'
  memory { params.genomeSize > 200000000 ? 50.GB : 16.GB }



  input:
  tuple val(sampleID), file(wig)


  output:
  tuple val(sampleID), file("*_FuzzSort.xls")

  script:
  """
  danpos.py dpos $wig -z 20 -e 1 \
   --width 10 --height 25 > $sampleID"_DANPOS_stats.txt"

   awk 'NR == 1; NR > 1 {print  \$0 | "sort -k6,6" }' result/pooled/*.xls > $sampleID"_FuzzSort.xls"

  """
}
