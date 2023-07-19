// load modules
//fastqc
include{fastqc; multiqc} from '../modules/raw_qc'
//trimming
include{trimgalore} from '../modules/trim'
//alignment to ref genome
include{alignment} from '../modules/align'
//qualimap after alignment
include{qualimap} from '../modules/qualimap'
//InsertSize_Histogram
include{InsertSize_Histogram} from '../modules/InsertSize_Histogram'
// filtering sizes using alignmentSieve
include{sieve} from '../modules/alignmentsieve'
// prepare for DANPOS
include{chrsize} from '../modules/prepareDANPOS'
// DANPOS run
include{danpos} from '../modules/DANPOS'
//FragmentStatistics
include{statistics_read; statistics_plot} from '../modules/fragment_statistics'
//deeptools TSS
include{TSS_profile;TSS_profile_plot} from '../modules/get_TSS_profile'
//TSS Profile monoNucs
include{make_TSS_plots} from '../modules/make_TSS_plots'
//Nucleosome repeat length
include{NRL}from '../modules/NRL' 



workflow profiler{
        take:
        sampleSingle_ch
        samplePair_ch

        main:
        //fastqc of raw data
        fastqc(sampleSingle_ch)
        //trimming
        trimgalore(samplePair_ch)
        //alignment against reference with bowtie2 and subsequent MAPQC filtering
        alignment(trimgalore.out[0])
        //QC of alignment
        qualimap(alignment.out[1])
        //Fragment sizes
        InsertSize_Histogram(qualimap.out[0].collect())
        // get mono-nucleosome fragments
        sieve(alignment.out[3])
        //get chrom_Sizes
        chrsize(sieve.out[1].map{name,bam -> file(bam)}.collect(sort:true).map { it[0] })
        //get nucleosome profiles
        danpos(sieve.out[1], chrsize.out)
        //FragmentStatistics
        statistics_read(sieve.out[0].join(fastqc.out[2]).join(trimgalore.out[4])
        .join(alignment.out[2]).join(qualimap.out[1]))
        statistics_plot(statistics_read.out[0].collect())

        //QualityCheck
        multiqc(fastqc.out[0].mix(trimgalore.out[1]).mix(trimgalore.out[2])
        .mix(alignment.out[0]).mix(qualimap.out[0]).collect())

        //TSS_Profile_mono
        if(params.TSS){
                TSS_profile(danpos.out[0].collect())
                TSS_profile_plot(TSS_profile.out)
                make_TSS_plots(TSS_profile_plot.out)
        }

        //NRL
        NRL(sieve.out[2].map{id, bam, idx -> id}.toList(), sieve.out[2].map{id, bam, idx -> bam}.toList(),sieve.out[2].map{id, bam, idx -> idx}.toList())
        
}
