process TSS_profile{
  cpus = (Math.round(Runtime.runtime.availableProcessors()*0.8))
  memory { 2.GB * task.cpus }

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
