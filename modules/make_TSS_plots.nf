process make_TSS_plots{
  container 'uschwartz/r_nucmacc:v3.1'

  publishDir "${params.outDir}/RUN/02_TSS_profile", mode: 'copy'

  input:
  file(input)

  output:
  file("*.pdf")

  script:
  """
  make_TSS_plots_monoNucs.R $input
  """
}
