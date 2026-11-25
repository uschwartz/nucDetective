process make_TSS_plots{
  container 'uschwartz/r_nucmacc:v3.1'
  tag "$ref_point"

  // Inspector: RUN/07_METAGENE_PLOTS/<TSS|TES>/, Profiler: RUN/02_TSS_profile/
  publishDir path: { params.analysis == 'inspector'
                       ? "${params.outDir}/RUN/07_METAGENE_PLOTS/${ref_point}"
                       : "${params.outDir}/RUN/02_${ref_point}_profile" }, mode: 'copy'

  input:
  tuple val(ref_point), file(input)

  output:
  file("*.pdf")

  script:
  """
  make_TSS_plots_monoNucs.R $input $ref_point
  """
}
