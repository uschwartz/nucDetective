#!/usr/bin/env nextflow

 /*
 ===============================================================================
                      nextflow based MNase-seq analysis pipeline
 ===============================================================================
Authors:
Uwe Schwartz <uwe.schwartz@ur.de>
 -------------------------------------------------------------------------------
 */

// settings and help message
include { settings    } from './modules/setting'
include { helpMessage } from './modules/help'
// generate profiles
include { profiler    } from './workflows/profiler'
// get dynamic nucleosomes
include { inspector   } from './workflows/inspector'


// Typed parameters. Since Nextflow 26 (strict syntax) values given on the command
// line arrive as Strings ("--normalize_profiles false" -> "false", which is truthy).
// Declaring the type here makes Nextflow convert them. Defaults stay in
// nextflow.config (needed there for config-time settings) and in conf/test_*.config.
// Not declared on purpose:
//  - path options that default to false (TSS, blacklist, nrl_regions, scratch_dir):
//    they are either false or a path, which a single type cannot express
//  - genomeSize: may exceed Integer range (Long is not accepted on the CLI in
//    Nextflow 26.10); nextflow.config casts it with `as long`, modules only print it
params {
    help:               Boolean
    test:               Boolean
    publishBam:         Boolean
    publishBamFlt:      Boolean
    OmitPublishWig:     Boolean
    normalize_profiles: Boolean
    detect_shift:       Boolean
    regularity_ref:     Boolean
    MAPQC:              Integer
    minLen:             Integer
    maxLen:             Integer
    peaks_used:         Integer
    avgNRL:             Integer
    nucref_min_support: Integer
}


workflow {

    // help message
    if (params.help) {
        helpMessage()
        exit(0)
    }

    // validate parameters
    if (!(params.container_engine in ['docker', 'singularity'])) {
        error "Invalid container engine: '${params.container_engine}'. Please use either 'docker' or 'singularity'."
    }
    if (!(params.analysis in ['profiler', 'inspector'])) {
        error "Invalid analysis: '${params.analysis}'. Please use either 'profiler' or 'inspector'."
    }

    // show settings
    settings()

    // read csv file; paths in the toyData sample sheets are relative to the pipeline directory
    def prefix  = params.test ? params.project : ''
    def rows_ch = channel.fromPath(params.csvInput).splitCsv(header: true)

    if (params.analysis == 'profiler') {
        // channel for alignment: [name, fwd, rev]
        def samplePair_ch = rows_ch.map { row ->
            tuple(row.Sample_Name, file(prefix + row.path_fwdReads), file(prefix + row.path_revReads))
        }
        // channel for fastqc: [name, fwd] and [name, rev]
        def sampleSingle_ch = samplePair_ch.flatMap { name, fwd, rev ->
            [ tuple(name, fwd), tuple(name, rev) ]
        }

        profiler(sampleSingle_ch, samplePair_ch)
    }

    if (params.analysis == 'inspector') {
        def samples_wig_ch = rows_ch.map { row ->
            tuple(row.Sample_Name, file(prefix + row.path_wig), row.normalize)
        }
        // the README and toyData use 'Condition'; older sample sheets use 'condition'
        def condition_ch = rows_ch.map { row ->
            tuple(row.Sample_Name, row.Condition ?: row.condition)
        }

        inspector(samples_wig_ch, condition_ch)
    }
}
