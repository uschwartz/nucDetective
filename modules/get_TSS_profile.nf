// Profiles around a reference point of the transcript annotation (GTF).
// ref_point: 'TSS' (transcription start site) or 'TES' (transcription end site)
process TSS_profile{
  // resources (cpus/memory): see process withName:'TSS_profile' in nextflow.config
  tag "$ref_point"

  input:
  file(bw)
  tuple val(ref_point), file(gtf)

  output:
  tuple val(ref_point), file("computeMatrix2plot_mono.txt.gz")


  script:
  """
  computeMatrix reference-point -S $bw \
   -R $gtf \
   --referencePoint $ref_point \
   -o "computeMatrix2plot_mono.txt.gz" \
   -b 1500 -a 1500 --smartLabels -p $task.cpus


  """
}

process TSS_profile_plot{
  tag "$ref_point"

  input:
  tuple val(ref_point), file(computeMatrix2plot_mono)

  output:
  tuple val(ref_point), file("values_Profile_mono.txt")

  script:
  """
  plotProfile -m $computeMatrix2plot_mono \
       -out 'DefaultHeatmap_mono.png' \
       --outFileNameData 'values_Profile_mono.txt'
  """
}
