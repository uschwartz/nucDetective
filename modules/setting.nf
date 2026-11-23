/*
* Show settings at the beginning of a run
*/

def settings() {
    def general = """
        nucDetective   P I P E L I N E   v${workflow.manifest.version}
        =============================
        Path variables used in analysis
        csvInput:             ${params.csvInput}
        outDir:               ${params.outDir}
        genomeIdx:            ${params.genomeIdx}
        chrSizes:             ${params.chrSizes}

        General options
        analysis:             ${params.analysis}
        test:                 ${params.test}
        container_engine:     ${params.container_engine}
        scratch_dir:          ${params.scratch_dir}
        TSS:                  ${params.TSS}
        genomeSize:           ${params.genomeSize}
        """.stripIndent()

    def profiler = """
        'profiler' options
        blacklist:            ${params.blacklist}
        MAPQC:                ${params.MAPQC}
        minLen:               ${params.minLen}
        maxLen:               ${params.maxLen}
        peaks_used:           ${params.peaks_used}
        nrl_regions:          ${params.nrl_regions}
        publishBam:           ${params.publishBam}
        publishBamFlt:        ${params.publishBamFlt}
        OmitPublishWig:       ${params.OmitPublishWig}
        """.stripIndent()

    def inspector = """
        'inspector' options
        normalize_profiles:   ${params.normalize_profiles}
        nucref_min_support:   ${params.nucref_min_support}
        detect_shift:         ${params.detect_shift}
        regularity_ref:       ${params.regularity_ref}
        avgNRL:               ${params.avgNRL}
        """.stripIndent()

    def resources = """
        Resource limits (per task)
        max_cpus:             ${params.max_cpus}
        max_memory:           ${params.max_memory}
        max_time:             ${params.max_time}
        """.stripIndent()

    log.info general + (params.analysis == 'inspector' ? inspector : profiler) + resources
}
