process regularity{
    container 'leoschmutterer/regularity:v1.1'
    cpus = 5
    memory = { 4.GB * task.cpus }
    publishDir "${params.outDir}/RUN/07_REGULARITY/HighFuzzRegions/", mode: 'copy', pattern: "HighFuzRegions*"
    publishDir "${params.outDir}/RUN/07_REGULARITY/RollingWindows/", mode: 'copy', pattern: "*RollingWindow.bw"
    publishDir "${params.outDir}/RUN/07_REGULARITY/ResultTable/", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Cutoff/", mode: 'copy', pattern: "*.pdf"

    input: 
    val(sampleID)
    file(bw)

    output:
    
    file("*.pdf")
    file("*.bed")
    file("*.bw")
    file("*.tsv")
    
    script:
  """
  regularity.R $sampleID $bw $params.chrSizes
  """
}
