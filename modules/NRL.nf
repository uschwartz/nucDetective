process NRL {
    container = 'leoschmutterer/nrl:v2.1'
    memory = 8.GB
    publishDir "${params.outDir}/RUN/03_NRL/", mode: 'copy', pattern: "NRLs.pdf"
    publishDir "${params.outDir}/QC/09_NRL/${id}/", mode: 'copy', pattern: "*monoNuc.pdf"
    publishDir "${params.outDir}/RUN/03_NRL/", mode: 'copy', pattern: "*.csv"
    
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
