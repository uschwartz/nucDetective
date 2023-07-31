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
