/*
* Show settings at the beginning
*/

def settings() {
        println ''
         log.info """\

                  nucMACC   P I P E L I N E
                  =============================
                  Path variables used in analysis
                  csvInput :            ${params.csvInput}
                  outDir   :            ${params.outDir}
                  genomeIdx:            ${params.genomeIdx}
                  chrSizes:             ${params.chrSizes}

                  General options
                  analysis:             ${params.analysis}
                  blacklist:            ${params.blacklist}
                  genomeSize:           ${params.genomeSize}
                  TSS:                  ${params.TSS}
                  publishBam:           ${params.publishBam}
                  publishBamFlt:        ${params.publishBamFlt}
                  OmitPublishWig:       ${params.OmitPublishWig}

                  """.stripIndent()

        println ''
}
