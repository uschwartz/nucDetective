/*
* Help message
*/

def helpMessage() {
    println ''
    log.info """
    nucMACC   P I P E L I N E
    =============================
    Usage:

    nextflow run uschwartz/nucMACC --csvInput 'path2csvFile' --outDir 'path2outDir' --genomeIdx 'path2bowtie2_idx' --genomeSize 'eff. genome size' --genome 'path2ref_genome_fasta'

    Mandatory arguments:
      --csvInput        [string] Path to comma-separated file containing information about the samples in the experiment (see ./toyData/input.csv as example) template provided in ./input_template.csv
      --analysis        [string] select the workflow to execute; available workflows 'profiler' or 'inspector'  (default:profiler)
      --genomeSize      [integer] Effective genome size, defined as the length of the mappable genome. Used for normalisation (default: 162367812 (dm3))

    Mandatory arguments profiler:
      --genomeIdx       [string] Path and prefix of bowtie2 index (minus .X.bt2)

    Mandatory arguments inspector:
      --chrSizes        [string] Path to tab-delim file with 2 columns 1)chr name and 2)chr size

    optional arguments:
      --outDir          [string] Name of output directory, which will be created (default: ~/nucMACC_test/)
      --blacklist       [string] A BED file containing regions that should be excluded from all nucleosome analysis (default: false)
      --TSS             [string] A Transcript annotation file in GTF format to obtain TSS profiles (default: false)
      --publishBam      [boolean]if set, aligned bam files will be written to outDir (default: false)
      --publishBamFlt   [boolean]if set, size selected bam files will be written to outDir (default: false)
      --OmitPublishWig  [boolean]if set, wig files will not be written to outDir (default: false)

<<<<<<< Updated upstream
=======
    inspector optional arguments:
    --regularity_ref    [boolean] Use reference map as reference positions for regularity (default: true)
    --normalize_profiles [boolean] Normalize nucleosome profile prior to analysis (default: true)
>>>>>>> Stashed changes
     """.stripIndent()
     println ''
}
