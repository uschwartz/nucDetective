process regularity{
    container 'leoschmutterer/regularity:v1.1'
    cpus = 4
    memory = { 15.GB * task.cpus }
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/", mode: 'copy', pattern: "Irregular*"
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/RollingWindows/", mode: 'copy', pattern: "*RollingWindow.bw"
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/", mode: 'copy', pattern: "*.png"
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/PSD/", mode: 'copy', pattern: "PSD_avgNRL_*.bw"


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
    cpus = 4
    memory = { 15.GB * task.cpus }
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/", mode: 'copy', pattern: "Irregular*"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/RollingWindows/", mode: 'copy', pattern: "*RollingWindow.bw"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/", mode: 'copy', pattern: "*.png"
    publishDir "${params.outDir}/RUN/07_REGULARITY/Reference_pos/PSD/", mode: 'copy', pattern: "PSD_avgNRL_*.bw"

    input: 
    val(sampleID)
    file(bw)
    file(reference_positions)

    output:
    
    file("*.png")
    file("*.bed")
    file("*.bw")
    file("*.tsv")
    
    script:
  """
  regularity.R $sampleID $bw $params.chrSizes $reference_positions
  """
}
