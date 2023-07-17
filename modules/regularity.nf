process regularity{
    container 'leoschmutterer/regularity:v1.1'
    cpus = 5
    memory = { 4.GB * task.cpus }
    publishDir "${params.outDir}/RUN/07_REGULARITY", mode: 'copy'

    input: 
    val(sampleID)
    file(bw)

    output:
    
    file("*.pdf")
    file("*.bed")
    file("*.bw")
    
    script:
  """
  regularity.R $sampleID $bw $params.chrSizes
  """
}
