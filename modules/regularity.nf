process regularity{
  
  if(params.container_engine == 'docker'){
    containerOptions "-v \$(dirname ${params.chrSizes}):\$(dirname ${params.chrSizes})"
  }
  if(params.container_engine == 'singularity'){
    containerOptions "-B \$(dirname ${params.chrSizes}):\$(dirname ${params.chrSizes})"
  }  

    container 'leoschmutterer/regularity:v1.1'
    // resources: see process withName:'regularity|regularity_reference' in nextflow.config
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/", mode: 'copy', pattern: "Irregular*"
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/RollingWindows/", mode: 'copy', pattern: "*RollingWindow.bw"
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/", mode: 'copy', pattern: "*.png"
    publishDir "${params.outDir}/RUN/07_REGULARITY/All_regions/PSD/", mode: 'copy', pattern: "PSD_avgNRL_*.bw"

    input: 
    val(sampleID)
    file(bw)

    output:
    
    file("*.png")
    file("*.bed")
    file("*.bw")
    file("*.tsv")
    
    script:
  """
  REGULARITY_WORKERS=$task.cpus regularity.R $sampleID $bw $params.chrSizes $params.avgNRL
  """
}

process regularity_reference{
  if(params.container_engine == 'docker'){
    containerOptions "-v \$(dirname ${params.chrSizes}):\$(dirname ${params.chrSizes})"
  }
  if(params.container_engine == 'singularity'){
    containerOptions "-B \$(dirname ${params.chrSizes}):\$(dirname ${params.chrSizes})"
  }  

    container 'leoschmutterer/regularity:v1.1'
    // resources: see process withName:'regularity|regularity_reference' in nextflow.config
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
  REGULARITY_WORKERS=$task.cpus regularity.R $sampleID $bw $params.chrSizes $reference_positions $params.avgNRL
  """
}
