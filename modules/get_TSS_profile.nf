process TSS_profile{
  // resources (cpus/memory): see process withName:'TSS_profile' in nextflow.config

  input:
  file(bw)
  file(tss)

  output:
  file "computeMatrix2plot_mono.txt.gz"


  script:
  """
  computeMatrix reference-point -S $bw \
   -R $tss \
   --referencePoint TSS \
   -o "computeMatrix2plot_mono.txt.gz" \
   -b 1500 -a 1500 --smartLabels -p $task.cpus


  """
}

process TSS_profile_plot{


  input:
  file(computeMatrix2plot_mono)

  output:
  file "values_Profile_mono.txt"

  script:
  """
  plotProfile -m $computeMatrix2plot_mono \
       -out 'DefaultHeatmap_mono.png' \
       --outFileNameData 'values_Profile_mono.txt'
  """
}
