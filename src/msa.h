#ifndef MSA_H
#define MSA_H

#include "types.h"
#include "fasta.h"
#include "features_profile.h"
#include "profile.h"

#include <iostream>
#include <vector>


class Sequences;


namespace msa {
  /// performs the full multiple sequence alignment, returns aligned sequences
  std::vector<std::string> run_msa(fasta::FastaData fasta_data, double gap_open_pen,
                                   double gap_ext_pen, double end_pen,
                                   int domain_score, int motif_score,
                                   int phosph_score, int codon_length);


  /// performs the first round of alignments, / all vs query seq (first
  /// calculates profile / based only on the query seq, then / aligns all
  /// sequences and calculates / identity of each sequence to the query seq.)
  
  
  /* UNCOMMENT WHEN FeaturesProfile struct is finished */
  // std::vector<std::string> PerformMSAfirstRound(fasta::FastaData fasta_data,
  //     const ProfileMap& profile, FeaturesProfile& output_features_profile,
  //     double penalty, double endPenalty, double extensionPenalty,
  //     int codon_length, IdentitiesList& identities);
  //    ###########

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
  // void RemoveGaps(fasta::Sequence& alignment_with_lowercase,
  //                 fasta::Sequence& alignment_without_lowercase,
  //                 std::vector<fasta::Sequence>& alignment);
  //    ###########
  ///
  /// takes a sequence and profiles, returns an
  /// alignment of the two, with gaps cut out
  ///
  //    ###########
  // void AlignPairwise(fasta::Sequence& al_without_lower,
  //                    fasta::Sequence& al_with_lower,
  //                    fasta::Sequence& seq2, const ProfileMap& prf,
  //                    FeaturesProfile& featPrf,
  //                    double penalty, double endPenalty, double extensionPenalty,
  //                    int codon_length);
  //    ###########

  ///
  /// calculates identity with the query sequence
  /// @param alignedSequence sequence aligned to the profile with the gaps cut
  /// out (its length is equal to the profile's length)
  ///
  //    ###########
  // double CalcIdentity(const fasta::Sequence& alignedSequence);
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
