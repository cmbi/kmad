#ifndef MSA_H
#define MSA_H

#include "types.h"
#include "seq_data.h"
#include "features_profile.h"
#include "profile.h"

#include <iostream>
#include <vector>


class Sequences;


namespace msa {
  /// performs the full multiple sequence alignment, returns aligned sequences
  std::vector<std::string> run_msa(const seq_data::SequenceData& sequence_data,
                                   const f_config::FeatureSettingsMap& f_set,
                                   double gap_open_pen,
                                   double gap_ext_pen, double end_pen,
                                   int domain_score, int motif_score,
                                   int phosph_score, int codon_length);


  /// performs the first round of alignments, / all vs query seq (first
  /// calculates profile / based only on the query seq, then / aligns all
  /// sequences and calculates / identity of each sequence to the query seq.)
  
  
  std::vector<double> set_identities(
      const seq_data::SequenceData& sequence_data,
      const ProfileMap& profile, FeaturesProfile& f_profile,
      double penalty, double endPenalty, double extensionPenalty,
      int codon_length);

  ///
  /// performs next round of MSA (good for all rounds except for the first one
  /// - you need a profile)
  ///
  
  
  /* UNCOMMENT WHEN FeaturesProfile struct is finished */
  // void PerformMSAnextRound(std::vector<std::string>& prevAlignment,
  //                          const ProfileMap& profile,
  //                          FeaturesProfile& output_features_profile,
  //                          double penalty, double endPenalty,
  //                          double extensionPenalty,
  //                          double identityCutoff,
  //                          int codon_length, IdentitiesList& identities,
  //                          int& prev_alignments);
  //    ###########

  ///
  /// takes pairwise alignment, removes
  /// characters from the 2nd sequence that match gaps from 1st seq and returns
  /// vector<string> of 2 elements, where the 1st one is 2nd sequence with cut
  /// out chars and 2nd one
  /// is 2nd sequence with cut out chars and lowercase chars
  /// before and after that
  ///
  //    ###########
  fasta::SequenceList remove_gaps(const fasta::SequenceList& alignment);
  //    ###########
  ///
  /// takes a sequence and profiles, returns an
  /// alignment of the two, with gaps cut out
  ///
  fasta::SequenceList align_pairwise(const fasta::Sequence& input_sequence, 
                                     const ProfileMap& profile, 
                                     const FeaturesProfile& f_profile,
                                     double gap_open_pen, double end_pen,
                                     double gap_ext_pen, int codon_length);

  ///
  /// calculates identity with the query sequence
  /// @param alignedSequence sequence aligned to the profile with the gaps cut
  /// out (its length is equal to the profile's length)
  ///
  //    ###########
  double calc_identity(const fasta::Sequence& aligned_sequence, 
                       const fasta::Sequence& query_sequence);
  //    ###########

  //    ###########
  // void add_feature_indexes(FeaturesProfile& fprf);
  //    ###########
  ///
  /// count alignments that will be performed in this round
  ///
  //    ###########
  // int CountAlignments(double identity_cutoff, IdentitiesList& identities);
  //    ###########
}

#endif /* MSA_H */
