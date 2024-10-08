#!/usr/bin/env nextflow

 /*
 ===============================================================================
                      nextflow based MNase-seq analysis pipeline
 ===============================================================================
Authors:
Uwe Schwartz <uwe.schwartz@ur.de>
 -------------------------------------------------------------------------------
 */

nextflow.enable.dsl = 2

 //                           show settings
 if (!params.help) {
         include{settings} from './modules/setting'
         settings()
 }

 //                       help message
 // Show help message
 if (params.help) {
     include{helpMessage} from './modules/help'
     helpMessage()
     exit 0
 }

//                      workflow
// read csv file

if(params.analysis=='profiler'){
        // forward reads
        Channel
            .fromPath(params.csvInput)
            .splitCsv(header:true)
            .map{ row -> tuple(row.Sample_Name,file(row.path_fwdReads))}
            .set{samples_fwd_ch}
        // reverse reads
        Channel
              .fromPath(params.csvInput)
              .splitCsv(header:true)
              .map{ row -> tuple(row.Sample_Name,file(row.path_revReads))}
              .set{samples_rev_ch}
        //Channel for fastqc
        samples_fwd_ch.mix(samples_rev_ch).set{sampleSingle_ch}
        //Channel for alignment
        samples_fwd_ch.join(samples_rev_ch).set{samplePair_ch}

        
}

if(params.analysis=='inspector'){
        Channel
            .fromPath(params.csvInput)
            .splitCsv(header:true)
            .map{ row -> tuple(row.Sample_Name,file(row.path_wig), row.normalize)}
            .set{samples_wig_ch}

        Channel
            .fromPath(params.csvInput)
            .splitCsv(header:true)
            .map{ row -> tuple(row.Sample_Name, row.condition)}
            .set{condition_ch}
}

// load workflows
// generate profiles
include{profiler} from './workflows/profiler'
// get dynamic nucleosomes
include{inspector} from './workflows/inspector'

workflow{

  if(params.analysis=='profiler'){
          profiler(sampleSingle_ch,samplePair_ch)
  }

  if(params.analysis=='inspector'){
          inspector(samples_wig_ch,condition_ch)
  }

}
