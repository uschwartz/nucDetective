// Distribution of dynamic nucleosome features (regularity changes, occupancy changes,
// shifts, fuzziness) along genes. Requires a gene annotation in GTF format (--TSS).
process feature_profile{
  // resources (cpus/memory): see process withName:'feature_profile' in nextflow.config
  publishDir "${params.outDir}/RUN/07_METAGENE_PLOTS/FEATURE_DISTRIBUTION/", mode: 'copy', pattern: "{values_featureDistribution.txt,DefaultProfile_featureDistribution.png}"

  input:
  tuple val(labels), file(beds)
  file(gtf)
  file(chrSizes)

  output:
  file "values_featureDistribution.txt"
  file "DefaultProfile_featureDistribution.png"

  script:
  def features = [labels, beds instanceof List ? beds : [beds]].transpose().collect { l, b -> "${l}:${b}" }.join(' ')
  """
  feature_distribution_matrix.sh $gtf $chrSizes $task.cpus $features
  """
}

process plot_feature_distribution{
  container 'uschwartz/r_nucmacc:v3.1'

  publishDir "${params.outDir}/RUN/07_METAGENE_PLOTS/FEATURE_DISTRIBUTION/", mode: 'copy'

  input:
  file(values)

  output:
  file("featureDistribution.pdf")
  file("relative_featureDistribution.txt")

  script:
  """
  plot_feature_distribution.R $values 750 750
  """
}
