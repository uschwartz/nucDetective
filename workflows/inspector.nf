// load modules
//quantile normalization
include{quantNorm;ref_bw} from '../modules/quantileNorm'
//nucleosome calling
include{danpos_nuc} from '../modules/DANPOS_finalnucs'
//nucleosome format conversion
include{nucs2bed} from '../modules/nucs2bed'
//get referenceMap
include{reference_map;reference_map_all } from '../modules/reference_map'
//deeptools TSS
include{TSS_profile;TSS_profile_plot} from '../modules/get_TSS_profile'
//TSS Profile monoNucs
include{make_TSS_plots} from '../modules/make_TSS_plots'



workflow inspector{
        take:
        wig_channel
        cond_channel

        main:
        //get samples for normalisation
        wig_channel.map{
           name,wig,norm -> tuple(name,wig,norm.toLowerCase())
         }.filter { it[2]!='reference'}.set{sample_ch}
         //get reference sample
         wig_channel.map{
            name,wig,norm -> tuple(name,wig,norm.toLowerCase())
          }.filter { it[2]=='reference'}.set{ref_ch}
          // quantile norm
          quantNorm(sample_ch,ref_ch)
          ref_bw(ref_ch)

          //nuc calling
          quantNorm.out[1].mix(ref_ch.map{name,wig,norm -> tuple(name,wig)}).set{norm_ch}
          danpos_nuc(norm_ch)
          //get top20 percent
          nucs2bed(danpos_nuc.out)
          //get referenceMap
          cond_channel.toSortedList{ a, b -> a[1] <=> b[1]}.transpose().first()
          .set{sorted_ch}
          reference_map(nucs2bed.out[0].collect(), sorted_ch)
          reference_map_all(nucs2bed.out[1].collect(), sorted_ch)

          //TSS_Profile
          if(params.TSS){
                  quantNorm.out[0].mix(ref_bw.out[0])
                  .map{name,bw -> bw}.set{bw_ch}
                  TSS_profile(bw_ch.collect())
                  TSS_profile_plot(TSS_profile.out)
                  make_TSS_plots(TSS_profile_plot.out)
          }


}
