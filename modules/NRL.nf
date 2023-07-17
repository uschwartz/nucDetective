process NRL {
    container = 'leoschmutterer/nrl:v2.1'
    memory = 8.GB
    publishDir "${params.outDir}/RUN/03_NRL/", mode: 'copy'

    input:
    val(id)
    file(bam)
    file(idx)

    output:
    file ("*.pdf")
    file (".csv")

   script:
   """
   NRL.R $id $bam
   """

}
