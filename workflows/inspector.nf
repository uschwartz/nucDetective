// load modules
//quantile normalization
include{quantNorm;ref_bw;wig_to_bw} from '../modules/quantileNorm'
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
//get peak scores
include{score_peaks;  pca; correlation} from '../modules/score_peaks'
//get dynamic nucleosomes
include{occupancy; shift} from '../modules/DynNucs'
//get fuzziness of nucleosomes
include{fuzziness} from '../modules/fuzziness'
//get regularity
include{regularity} from '../modules/regularity'


workflow inspector{
        take:
        wig_channel
        cond_channel

        main:
        //get samples for normalisation
        if(params.normalize_profiles){
        wig_channel.map{
           name,wig,norm -> tuple(name,wig,norm.toLowerCase())
         }.filter { it[2]!='reference'}.set{sample_ch}
         //get reference sample
         wig_channel.map{
            name,wig,norm -> tuple(name,wig,norm.toLowerCase())
          }.filter { it[2]=='reference'}.set{ref_ch}
          sample_ch.combine(ref_ch)
          // quantile norm
          quantNorm(sample_ch.combine(ref_ch))
          ref_bw(ref_ch)

          //nuc calling
          quantNorm.out[1].mix(ref_ch.map{name,wig,norm -> tuple(name,wig)}).set{norm_ch}
          danpos_nuc(norm_ch)

          quantNorm.out[0].mix(ref_bw.out[0]).map{name,bw -> bw}.set{bw_ch}
          } 
          else {
          //get samples
          wig_channel.map{name,wig,norm -> tuple(name,wig)}.set{sample_ch}
          //nuc calling
          danpos_nuc(sample_ch)
          //get bw
          wig_to_bw(sample_ch)
          wig_to_bw.out.map{name,bw -> bw}.set{bw_ch}
          }

          //get top20 percent
          nucs2bed(danpos_nuc.out)
          //get referenceMap
          cond_channel.toSortedList{ a, b -> a[1] <=> b[1]}.transpose().first()
          .set{sorted_ch}
          reference_map(nucs2bed.out[0].collect(), sorted_ch)
          reference_map_all(nucs2bed.out[1].collect(), sorted_ch)
         
          
          

          //TSS_Profile
          if(params.TSS){
<<<<<<< Updated upstream
                  TSS_profile(bw_ch.collect())
=======
               tss_ch = Channel.fromPath(params.TSS)
        
                  TSS_profile(bw_ch.collect(),tss_ch)
>>>>>>> Stashed changes
                  TSS_profile_plot(TSS_profile.out)
                  make_TSS_plots(TSS_profile_plot.out)
          }

          //get score under nucleosomes positions
          score_peaks(bw_ch.collect(), reference_map.out)
          //exploratory data analysis
          pca(score_peaks.out[1])
          correlation(score_peaks.out[1])

          // dynamic nucleosomes
          occupancy(score_peaks.out[0])
          shift(reference_map.out, bw_ch.collect())

          //fuzziness of nucleosomes
          fuzziness(reference_map.out, nucs2bed.out[0].collect())

<<<<<<< Updated upstream
          //regularity of nucleosomes  
          regularity(quantNorm.out[0].mix(ref_bw.out[0]).map{id, bw -> id}.toSortedList(),
          quantNorm.out[0].mix(ref_bw.out[0]).map{id, bw -> bw}.toSortedList())
          

=======
          //regularity of nucleosomes
          if (params.normalize_profiles){
          if(params.regularity_ref){
                regularity_reference(quantNorm.out[0].mix(ref_bw.out[0]).map{id, bw -> id}.toSortedList(), quantNorm.out[0].mix(ref_bw.out[0]).map{id, bw -> bw}.toSortedList(), reference_map.out)
          } else {
                 regularity(quantNorm.out[0].mix(ref_bw.out[0]).map{id, bw -> id}.toSortedList(), quantNorm.out[0].mix(ref_bw.out[0]).map{id, bw -> bw}.toSortedList())       
          }
          }
          else{
            if(params.regularity_ref){
                regularity_reference(wig_to_bw.out.map{id, bw -> id}.toSortedList(), wig_to_bw.out.map{id, bw -> bw}.toSortedList(), reference_map.out)
          } else {
                 regularity(wig_to_bw.out.map{id, bw -> id}.toSortedList(), wig_to_bw.out.map{id, bw -> bw}.toSortedList())       
          }
          }
>>>>>>> Stashed changes
}
