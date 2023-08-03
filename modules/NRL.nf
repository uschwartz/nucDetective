process NRL {
    container = 'leoschmutterer/nrl:v2.1'
    memory = 4.GB
    publishDir "${params.outDir}/QC/09_NRL/", mode: 'copy', pattern: "*monoNuc.pdf"
    
    input:
    val(id)
    file(bam)
    file(idx)

    output:
    file ("*.pdf")
    file ("*.csv")

   script:
   """
   NRL.R $id $bam $params.peaks_used
   """

}

process NRL_overview {
    container = 'leoschmutterer/nrl:v2.1'
    memory = 4.GB
    publishDir "${params.outDir}/RUN/03_NRL/", mode: 'copy'

    input:
    file (csv)

    output:
    file ("*.pdf")
    file ("*.csv")

    script:
    """
    Compare_NRL.R
    """
}

process NRL_regions {
	container = 'leoschmutterer/nrl_regions:v1.1'
	label 'big'
	memory '16.GB'
	publishDir "${params.outDir}/RUN/03_NRL/NRL_Regions/Single/", mode: 'copy', pattern: "*.bed"
    publishDir "${params.outDir}/RUN/03_NRL/NRL_Regions/Single/", mode: 'copy', pattern: "*.pdf"
	
    input:
    val(id)
    file(bam)
    file(idx)

    output:
    file ("*.pdf")
    file ("*.csv")
    file ("*.bed")

   script:
   """
   NRL_regions.R $id $bam $params.nrl_regions $params.peaks_used
   """

}

process compare_NRL_regions {
    container = 'leoschmutterer/nrl_regions:v1.1'
    publishDir "${params.outDir}/RUN/03_NRL/NRL_Regions/Summary/", mode: 'copy'

    input:
    file(values)

    output:
    file ("*.csv")
    file("*.pdf")

    script:
   """
   NRL_region_comparison.R 
   """
}
