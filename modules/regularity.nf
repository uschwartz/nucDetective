process regularity{
    container 'leoschmutterer/regularity:v1.1'
    cpus = 5
    memory = { 4.GB * task.cpus }
    publishDir "${params.outDir}/RUN/07_REGULARITY/Irregular_Regions/", mode: 'copy', pattern: "Irregular*"
    publishDir "${params.outDir}/RUN/07_REGULARITY/RollingWindows/", mode: 'copy', pattern: "*RollingWindow.bw"
    publishDir "${params.outDir}/RUN/07_REGULARITY/ResultTable/", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Cutoff/", mode: 'copy', pattern: "*.pdf"
    publishDir "${params.outDir}/RUN/07_REGULARITY/PSD/", mode: 'copy', pattern: "PSD_avgNRL_*.bw"

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

process regularity_reference{
    container 'leoschmutterer/regularity:v1.1'
    cpus = 5
    memory = { 4.GB * task.cpus }
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/Irregular_Regions/", mode: 'copy', pattern: "Irregular*"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/RollingWindows/", mode: 'copy', pattern: "*RollingWindow.bw"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/ResultTable/", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/Cutoff/", mode: 'copy', pattern: "*.pdf"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/PSD/", mode: 'copy', pattern: "PSD_avgNRL_*.bw"

    input: 
    val(sampleID)
    file(bw)
    file(reference_positions)

    output:
    
    file("*.pdf")
    file("*.bed")
    file("*.bw")
    file("*.tsv")
    
    script:
  """
  regularity.R $sampleID $bw $params.chrSizes $reference_positions
  """
}
