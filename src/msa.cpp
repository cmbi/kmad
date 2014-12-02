#include "msa.h"
#include "txtproc.h"
#include "features_profile.h"
#include "profile.h"
#include "sequences.h"

#include <iostream>
#include <vector>
#include <tuple>


StringSequences msa::run_msa(Sequences sequences,
                              std::string conf_filename,
                              double gap_open_pen,
                              double gap_ext_pen,
                              double end_pen,
                              double lcr_mod,
                              int domain_score, 
                              int motif_score,
                              int phosph_score,
                              int codon_length,
                              IDsList motifs_ids,
                              ProbsList motifs_probs){

      Profile prf;
      FeaturesProfile fprf(domain_score, phosph_score, motif_score, lcr_mod, 
                           motifs_ids, motifs_probs);
      if (!conf_filename.empty()){
        txtproc::process_conf_file(conf_filename, fprf, sequences);
      }
      std::vector<double> identities;
      //first round of the alignment - all vs 1st
      std::vector<std::string> multipleAlignment;
      multipleAlignment = sequences.PerformMSAfirstRound(prf, fprf, 
                                                         gap_open_pen, 
                                                         end_pen, 
                                                         gap_ext_pen, 
                                                         codon_length, 
                                                         identities);
      StringSequences alignment;
      int prev_alignments = 0;
      for (int i = 8; i >= 0; i--){
        double cutoff = double(i)/10;
        sequences.PerformMSAnextRound(alignment, prf, fprf, gap_open_pen, 
                                      end_pen, gap_ext_pen,
                                      cutoff, codon_length, 
                                      identities, prev_alignments);
        //prev_alignments - number of alignments performed in the previous 
        //rounds - to omit this round if the number of aligned sequences is the
        //same as in the previous round
      }
      prev_alignments = 0;  // to align (again) all sequences to the profile
      sequences.PerformMSAnextRound(alignment, prf, fprf, 
                                    gap_open_pen, end_pen, gap_ext_pen, 0, 
                                    codon_length, identities, 
                                    prev_alignments);
      return alignment;
}
