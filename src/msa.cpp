#include "msa.h"

#include "features_profile.h"
#include "profile.h"
#include "scoring_matrix.h"
#include "substitution_matrix.h"
#include "vec_util.h"

#include <boost/algorithm/string.hpp>    
#include <boost/filesystem.hpp>

#include <iostream>
#include <vector>
#include <tuple>


fasta::SequenceList msa::run_msa(const seq_data::SequenceData& sequence_data,
                                 const f_config::FeatureSettingsMap& f_set,
                                 double gap_open_pen, double gap_ext_pen,
                                 double end_pen, int domain_modifier,
                                 int motif_modifier, int ptm_modifier,
                                 int codon_length)
{
      FeaturesProfile f_profile(sequence_data.feature_list, domain_modifier,
                                ptm_modifier, motif_modifier,
                                sequence_data.probabilities);
      // query_seq_list - the profile are built only based on the first
      // sequence
      fasta::SequenceList query_seq_list = {sequence_data.sequences[0]};
      ProfileMap profile = create_profile(query_seq_list);
      f_profile.create_score_features_profile(query_seq_list, f_set);

      // first round of the alignment - all vs 1st
      // TODO: first round is used to set identities. Return value is ignored.
      // Can this be improved?
      std::vector<double> identities;
      identities = msa::set_identities(sequence_data, profile, f_profile,
                                       gap_open_pen, end_pen, 
                                       gap_ext_pen, codon_length);


      fasta::SequenceList alignment;
      int prev_alignments = 0;
      double cutoff = 0;
      for (int i = 8; i >= 0; i--) {
        cutoff = double(i) / 10;
        alignment = msa::perform_msa_round(profile, f_profile, gap_open_pen,
                                           end_pen, gap_ext_pen, cutoff, 
                                           codon_length, identities, 
                                           prev_alignments);
        //prev_alignments - number of alignments performed in the previous
        //rounds - to omit this round if the number of aligned sequences is the
        //same as in the previous round
      }
      prev_alignments = 0;  // to align (again) all sequences to the profile
      cutoff = 0;
      alignment = msa::perform_msa_round(profile, f_profile, gap_open_pen, 
                                         end_pen, gap_ext_pen, cutoff,
                                         codon_length, identities,
                                         prev_alignments);
      return alignment;
}


std::vector<double> msa::set_identities(
    const seq_data::SequenceData& sequence_data, const ProfileMap& profile,
    FeaturesProfile& f_profile, double gap_open_pen,
    double end_pen, double gap_ext_pen, int codon_length)
{
  // identity of the 1st one to itself
  // to build the first profile based only on the first seqeunce
  std::vector<double> identities = {1.};

  //pairwise alignment without lowercase characters
  fasta::Sequence aligned_seq_uppercase; 
  //pairwise alignment with lowercase characters where chars were removed
  fasta::Sequence aligned_seq_with_lower; 

  for (size_t i = 1; i < sequence_data.sequences.size(); ++i) {
    // aligned_sequence: vector, first element is the aligned sequence only
    // with uppercase characters, second one is with lowercase where the gaps
    // were cut out
    fasta::SequenceList aligned_sequence = msa::align_pairwise(
        sequence_data.sequences[i], profile,
        f_profile, gap_open_pen, end_pen, 
        gap_ext_pen, codon_length);

    double identity = msa::calc_identity(aligned_sequence[0],
                                         sequence_data.sequences[0]);
    identities.push_back(identity);
  }
  return identities;
}


double msa::calc_identity(const fasta::Sequence& aligned_sequence,
                          const fasta::Sequence& query_sequence) {
  double identical_residues = 0;
  assert(aligned_sequence.residues.size() == query_sequence.residues.size());
  for (unsigned i = 0; i < aligned_sequence.residues.size(); i++) {
    if (aligned_sequence.residues[i].codon[0] 
        == query_sequence.residues[i].codon[0]) {
      identical_residues++;
    }
  }
  return identical_residues / double(query_sequence.residues.size());
}


fasta::SequenceList msa::remove_gaps(const fasta::SequenceList& alignment) {
  fasta::Sequence new_seq;
  fasta::SequenceList aligned_seq = {new_seq, new_seq};
  char gap = '-';
  bool lower_flag = false;
  for (size_t i = 0; i < alignment[0].residues.size(); i++) {
    if (alignment[0].residues[i].codon[0] == gap) {
      if (aligned_seq[1].residues.size() > 0) {
        //change previous character to lowercase
        int last_index = aligned_seq[1].residues.size() - 1;
        aligned_seq[1].residues[last_index].codon[0] = tolower(
            aligned_seq[1].residues[last_index].codon[0]);
      }
      // flag to true so that the next character is also lowercase
      lower_flag = true;
    } else {
      if (lower_flag) {   //lowercase char
        // add lowercase char to the alignment with lowercases
        fasta::Residue new_residue = alignment[1].residues[i];
        new_residue.codon[0] = tolower(new_residue.codon[0]);
        aligned_seq[1].residues.push_back(new_residue);
        // add uppercase alignment to the alignment without lowercases
        aligned_seq[0].residues.push_back(alignment[1].residues[i]);
        lower_flag = false;
      } else {
        // uppercase char
        // adds the same uppercase char to both alignments (with lowercases and
        // without lowercases)
        aligned_seq[0].residues.push_back(alignment[1].residues[i]);
        aligned_seq[1].residues.push_back(alignment[1].residues[i]);
      }
    }
  }
  return aligned_seq;
}


fasta::SequenceList msa::align_pairwise(const fasta::Sequence& input_sequence,
                                        const ProfileMap& profile,
                                        const FeaturesProfile& f_profile,
                                        double gap_open_pen, double end_pen,
                                        double gap_ext_pen,
                                        int codon_length) {
  // int profile_length = profile.begin()->second.size();
  fasta::SequenceList alignment;
  //ScoringMatrix scores(profile_length, seq2.size(), gap_open_pen,
                       //end_pen, gap_ext_pen);
  //scores.CalculateScores(seq2, profile, feat_prf,
                         //codon_length);
  //scores.PerformNWAlignment(alignment, seq2, profile, feat_prf,
                            //codon_length);
  fasta::SequenceList aligned_sequence;
  aligned_sequence = remove_gaps(alignment);
  return aligned_sequence;
}


// TODO: Implement
fasta::SequenceList msa::perform_msa_round(
    const ProfileMap& profile,
    const FeaturesProfile& f_profile,
    double gap_open_pen,
    double end_pen,
    double gap_ext_pen,
    double identity_cutoff,
    int codon_length,
    IdentitiesList& identities,
    int& prev_alignments)
{
  fasta::SequenceList result;
  //int next_alignments = CountAlignments(identity_cutoff, identities);
  //if (next_alignments > prev_alignments) {
    ////working alignment - without lowercase around cut out residues
    ////would make latter aligning more complicated
    //std::vector<fasta::Sequence> alignment_without_lowercase;
    ////lowercase before and after cut out residues -- final result
    //std::vector<fasta::Sequence> alignment_with_lowercase;
    //alignment_without_lowercase.push_back(m_sequences_aa[0]);
    //alignment_with_lowercase.push_back(m_sequences_aa[0]);
    //// tmp pairwise alignment (and so is al_with_lower)
    //fasta::Sequence al_without_lower;
    //fasta::Sequence al_with_lower;
    //for (int i = 1; i < m_seq_nr; i++) {
      //if (identities[i] > identity_cutoff) {
        //// NW alignment of the ith seq against the profile
        //msa::AlignPairwise(al_without_lower, al_with_lower, m_sequences_aa[i],
                      //profile, output_features_profile, gap_open_pen,
                      //end_pen, gap_ext_pen, codon_length);
        //alignment_without_lowercase.push_back(al_without_lower);
        //alignment_with_lowercase.push_back(al_with_lower);
      //}
    //}
    ////create features profile based on the 1st seq
    //output_features_profile.CreateProfile(alignment_without_lowercase,
                                          //codon_length);
    //profile.ProcessProfile(alignment_without_lowercase);
    //prev_alignment = vec_util::Flatten(alignment_with_lowercase);
    ////update number of performed alignments
    //prev_alignments = next_alignments;
  //}
  return result;
}


// int msa::CountAlignments(double identity_cutoff,
//                          IdentitiesList& identities) {
//   int count = 0;
//   for (auto& item: identities) {
//     if (item > identity_cutoff) {
//       count++;
//     }
//   }
//   return count;
// }
