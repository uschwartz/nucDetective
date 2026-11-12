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


workflow {

    // help message
    if (params.help) {
        helpMessage()
        return
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
