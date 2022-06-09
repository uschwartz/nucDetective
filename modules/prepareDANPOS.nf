process chrsize{

  input:
  file(bam)
  output:
  file("chrom_Sizes.txt")

  script:
  """
  samtools view -H $bam \
  | awk -v OFS='\t' '/^@SQ/ {print \$2,\$3}' \
  | sed  's/.N://g' >chrom_Sizes.txt
  """
}
