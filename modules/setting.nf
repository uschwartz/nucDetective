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


                  'profiler' options
                  
                  minLen:              ${params.minLen}
                  maxLen:              ${params.maxLen}
                  MAPQC:                ${params.MAPQC}
                  peaks_used:           ${params.peaks_used}
		  nrl_regions:		${params.nrl_regions}
		  
                  'inspector' options:
                  normalize_profiles:   ${params.normalize_profiles}
                  regularity_ref:       ${params.regularity_ref}

                  """.stripIndent()

        println ''
}
