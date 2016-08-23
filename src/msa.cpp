#include "msa.h"

#include "feature_scores.h"
#include "optimizer.h"
#include "profile.h"
#include "scoring_matrix.h"

#include <boost/filesystem.hpp>
#include <vector>
#include <iostream>


std::vector<fasta::SequenceList> msa::run_msa(
    const fasta::FastaData& fasta_data,
    const f_config::FeatureSettingsMap& f_set,
    double gap_open_pen, double gap_ext_pen,
    double end_pen, double domain_modifier,
    double motif_modifier, double ptm_modifier,
    double strct_modifier,
    int codon_length, bool one_round,
    const std::string& sbst_mat, const bool first_gapped, const bool optimize,
    const bool fade_out, const bool no_feat)
{
      FeatureScores f_profile(fasta_data.feature_list, domain_modifier,
                              ptm_modifier, motif_modifier, strct_modifier,
                              fasta_data.probabilities);
      // query_seq_list - the profiles are built based only on the first
      // sequence
      fasta::SequenceList query_seq_list = {fasta_data.sequences[0]};
      auto profile = profile::create_score_profile(query_seq_list, sbst_mat);
      std::vector<double> identities = {1};
      if (!no_feat) {
        f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
      }

      // Align all sequences vs first to determine the identities
      identities = msa::set_identities(fasta_data, profile, f_profile,
                                       gap_open_pen, end_pen,
                                       gap_ext_pen, codon_length, no_feat);


      std::vector<fasta::SequenceList> alignment;
      auto alignments_number = 0;
      double cutoff = 0;
      // pointer to the function performing single round of msa, can be either
      // for gapped or ungapped first sequence
      std::vector<fasta::SequenceList> (*perform_msa_round_ptr)(
        const fasta::FastaData& fasta_data,
        const fasta::FastaData& fasta_data_alignment,
        const profile::ProfileMap& profile,
        const FeatureScores& f_profile,
        double gap_open_pen, double end_pen,
        double gap_ext_pen,
        double identity_cutoff,
        int codon_length,
        const std::vector<double>& identities,
        int& prev_alignments,
        const f_config::FeatureSettingsMap& f_set,
        std::vector<fasta::SequenceList> previous_alignment, int refine_seq,
        const bool no_feat);
      if (first_gapped) {
        perform_msa_round_ptr = msa::perform_msa_round_gapped;
      } else {
        perform_msa_round_ptr = msa::perform_msa_round_ungapped;
      }
      if (!one_round) {
        for (int i = 8; i >= 0; --i) {
          cutoff = double(i) / 10;
          auto prev_alignments = alignments_number;
          alignment = perform_msa_round_ptr(fasta_data, fasta_data, profile,
                                            f_profile, gap_open_pen,
                                            end_pen, gap_ext_pen, cutoff,
                                            codon_length, identities,
                                            alignments_number, f_set,
                                            alignment, 0, no_feat);
          // prev_alignments - number of alignments performed in the previous
          // round, needed not to update the profiles if number of aligned
          // sequences hasn't changed
          if (prev_alignments < alignments_number){
            if (!no_feat){
              f_profile.update_scores(alignment[0], f_set, identities,
                                      fade_out);
            }
            profile = profile::create_score_profile(alignment[0], sbst_mat);
          }
        }
      }
      // set alignments number to 0 to align (again)
      // all sequences to the profile
      auto iterations = 1;
      if (one_round) {
        iterations = 1;
      }
      for (int i = 0; i < iterations; ++i) {
        alignments_number = 0;
        cutoff = 0;
        alignment = perform_msa_round_ptr(fasta_data, fasta_data, profile,
                                          f_profile, gap_open_pen,
                                          end_pen, gap_ext_pen, cutoff,
                                          codon_length, identities,
                                          alignments_number, f_set,
                                          alignment, 0, no_feat);
        if (!no_feat) {
          f_profile.update_scores(alignment[0], f_set, identities, fade_out);
        }
        profile = profile::create_score_profile(alignment[0], sbst_mat);
      }
      if (optimize) {
        auto counter = 0;
        std::vector<fasta::SequenceList> previous;
        while (!msa::compare_alignments(previous, alignment)
                && counter < 15) {
          previous = alignment;
          alignment = optimizer::optimize_alignment(alignment, domain_modifier,
              motif_modifier, ptm_modifier, sbst_mat);
          ++counter;
        }
        if (!no_feat) {
          f_profile.update_scores(alignment[0], f_set, identities, fade_out);
        }
        profile = profile::create_score_profile(alignment[0], sbst_mat);
        alignments_number = 0;
        cutoff = 0;
        alignment = perform_msa_round_ptr(fasta_data,
                                          fasta_data, profile,
                                          f_profile, gap_open_pen,
                                          end_pen, gap_ext_pen, cutoff,
                                          codon_length, identities,
                                          alignments_number, f_set,
                                          alignment, 0, no_feat);
      }
      return alignment;
}


std::vector<fasta::SequenceList> msa::refine_alignment(
    const fasta::FastaData& fasta_data_plain,
    const fasta::FastaData& fasta_data_alignment,
    const f_config::FeatureSettingsMap& f_set,
    double gap_open_pen, double gap_ext_pen,
    double end_pen, double domain_modifier,
    double motif_modifier, double ptm_modifier,
    double strct_modifier,
    int codon_length, bool one_round,
    const std::string& sbst_mat, const bool first_gapped, const bool optimize,
    const bool fade_out, int refine_seq, const bool no_feat)
{
      fasta::SequenceList query_seq = {fasta_data_plain.sequences[0]};
      auto profile_single = profile::create_score_profile(query_seq, sbst_mat);
      FeatureScores f_profile_single(
              fasta_data_plain.feature_list, domain_modifier, ptm_modifier,
              motif_modifier, strct_modifier, fasta_data_plain.probabilities
      );
      std::vector<double> identities = {1};
      f_profile_single.update_scores(query_seq, f_set, identities, fade_out);
      identities = msa::set_identities(fasta_data_plain, profile_single,
                                       f_profile_single, gap_open_pen, end_pen,
                                       gap_ext_pen, codon_length, no_feat);
      FeatureScores f_profile(fasta_data_alignment.feature_list,
                              domain_modifier, ptm_modifier, motif_modifier,
                              strct_modifier,
                              fasta_data_alignment.probabilities);
      // query_seq_list - the profiles are built based on all sequences

      fasta::SequenceList query_seq_list = {};
      for (int i = 0; i < refine_seq; ++i) {
        query_seq_list.push_back(fasta_data_alignment.sequences[i]);
      }
      // fasta::SequenceList query_seq_list = fasta_data_alignment.sequences;
      auto profile = profile::create_score_profile(query_seq_list, sbst_mat);
      if (!no_feat) {
        f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
      }

      // Align all sequences vs first to determine the identities

      std::vector<fasta::SequenceList> alignment;
      auto alignments_number = 0;
      double cutoff = 0;
      // pointer to the function performing single round of msa, can be either
      // for gapped or ungapped first sequence
      std::vector<fasta::SequenceList> (*perform_msa_round_ptr)(
        const fasta::FastaData& fasta_data,
        const fasta::FastaData& fasta_data_alignment,
        const profile::ProfileMap& profile,
        const FeatureScores& f_profile,
        double gap_open_pen, double end_pen,
        double gap_ext_pen,
        double identity_cutoff,
        int codon_length,
        const std::vector<double>& identities,
        int& prev_alignments,
        const f_config::FeatureSettingsMap& f_set,
        std::vector<fasta::SequenceList> previous_alignment,
        int refine_seq, const bool no_feat);
      if (first_gapped) {
        perform_msa_round_ptr = msa::perform_msa_round_gapped;
      } else {
        perform_msa_round_ptr = msa::perform_msa_round_ungapped;
      }
      // set alignments number to 0 to align (again)
      // all sequences to the profile
      alignments_number = 0;
      cutoff = 0;
      alignment = perform_msa_round_ptr(fasta_data_plain,
                                        fasta_data_alignment, profile,
                                        f_profile, gap_open_pen,
                                        end_pen, gap_ext_pen, cutoff,
                                        codon_length, identities,
                                        alignments_number, f_set, alignment,
                                        refine_seq, no_feat);
      if (!no_feat) {
        f_profile.update_scores(alignment[0], f_set, identities, fade_out);
      }
      profile = profile::create_score_profile(alignment[0], sbst_mat);
      // set alignments number to 0 to align (again)
      // all sequences to the profile
      alignments_number = 0;
      cutoff = 0;
      alignment = perform_msa_round_ptr(fasta_data_plain,
                                        fasta_data_alignment, profile,
                                        f_profile, gap_open_pen,
                                        end_pen, gap_ext_pen, cutoff,
                                        codon_length, identities,
                                        alignments_number, f_set,
                                        alignment, refine_seq, no_feat);
      if (optimize) {
        auto counter = 0;
        std::vector<fasta::SequenceList> previous;
        while (!compare_alignments(previous, alignment) && counter < 15) {
          previous = alignment;
          alignment = optimizer::optimize_alignment(alignment, domain_modifier,
              motif_modifier, ptm_modifier, sbst_mat);
          ++counter;
        }
      }
      return alignment;
}


std::vector<double> msa::set_identities(
    const fasta::FastaData& fasta_data,
    const profile::ProfileMap& profile,
    FeatureScores& f_profile, double gap_open_pen,
    double end_pen, double gap_ext_pen, int codon_length,
    const bool no_feat)
{
  // identity of the 1st one to itself
  std::vector<double> identities = {1.};
  //pairwise alignment without lowercase characters
  fasta::Sequence aligned_seq_uppercase;
  //pairwise alignment with lowercase characters where chars were removed
  fasta::Sequence aligned_seq_with_lower;
  auto gapped = true;
  for (size_t i = 1; i < fasta_data.sequences.size(); ++i) {
    // aligned_sequence: vector
    // first element is a dummy polyA sequence to indicate where are the gaps
    // in the profile,
    // second one is with lowercase where the gaps
    // were cut out
    fasta::SequenceList aligned_sequence = msa::align_pairwise(
        fasta_data.sequences[i], profile,
        f_profile, gap_open_pen, end_pen,
        gap_ext_pen, codon_length, gapped, no_feat);


    auto identity = msa::calc_identity(
            aligned_sequence[0], aligned_sequence[1], fasta_data.sequences[0]);
    identities.push_back(identity);
  }
  return identities;
}

double msa::calc_identity(const fasta::Sequence& dummy_sequence,
                          const fasta::Sequence& aligned_sequence,
                          const fasta::Sequence& query_sequence) {
  double identical_residues = 0;
  auto gap_count = 0;
  // sequences should be aligned (therefore lengths should be equal)
  assert(aligned_sequence.residues.size() == dummy_sequence.residues.size());
  for (unsigned i = 0; i < aligned_sequence.residues.size(); ++i) {
    if (dummy_sequence.residues[i].codon[0] == '-') {
      ++gap_count;
    } else if (aligned_sequence.residues[i].codon[0]
        == query_sequence.residues[i - gap_count].codon[0]) {
      ++identical_residues;
    }
  }
  return identical_residues / double(dummy_sequence.residues.size());
}


fasta::SequenceList msa::remove_gaps(const fasta::SequenceList& alignment) {
  fasta::Sequence new_seq;
  fasta::SequenceList aligned_seq = {new_seq, new_seq};
  auto gap = '-';
  auto lower_flag = false;
  for (size_t i = 0; i < alignment[0].residues.size(); ++i) {
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
                                        const profile::ProfileMap& profile,
                                        const FeatureScores& f_profile,
                                        double gap_open_pen, double end_pen,
                                        double gap_ext_pen,
                                        int codon_length,
                                        const bool first_gapped,
                                        const bool no_feat) {
  auto profile_length = profile.begin()->second.size();
  ScoringMatrix scores(profile_length, input_sequence.residues.size(),
                       gap_open_pen, end_pen, gap_ext_pen, no_feat);
  scores.calculate_scores(input_sequence, profile, f_profile, codon_length);
  fasta::SequenceList alignment;
  alignment = scores.backtrace_alignment_path(input_sequence,
                                              profile, f_profile,
                                              codon_length);
  if (!first_gapped) {
    alignment = msa::remove_gaps(alignment);
  }
  return alignment;
}


std::vector<fasta::SequenceList> msa::perform_msa_round_ungapped(
    const fasta::FastaData& fasta_data,
    const fasta::FastaData& fasta_data_alignment,
    const profile::ProfileMap& profile,
    const FeatureScores& f_profile,
    double gap_open_pen,
    double end_pen,
    double gap_ext_pen,
    double identity_cutoff,
    int codon_length,
    const std::vector<double>& identities,
    int& prev_alignments,
    const f_config::FeatureSettingsMap& f_set,
    std::vector<fasta::SequenceList> previous_alignment,
    int refine_seq, const bool no_feat)
{
  std::vector<fasta::SequenceList> alignment = {{fasta_data.sequences[0]},
                                                {fasta_data.sequences[0]}};
  size_t start = 1;
  if (refine_seq != 0) {
    start = refine_seq;
    for (int i = 1; i < refine_seq; ++i) {
      alignment[0].push_back(fasta_data_alignment.sequences[i]);
      alignment[1].push_back(fasta_data_alignment.sequences[i]);
    }
  }
  auto next_alignments = count_alignments(identity_cutoff, identities);
  if (next_alignments > prev_alignments) {
    fasta::SequenceList aligned_seq;
    for (size_t i = start; i < fasta_data.sequences.size(); ++i) {
      if (identities[i] >= identity_cutoff) {
        // NW alignment of the ith seq against the profile
        aligned_seq = msa::align_pairwise(fasta_data.sequences[i],
                                          profile, f_profile, gap_open_pen,
                                          end_pen, gap_ext_pen, codon_length,
                                          false, no_feat);
        // if not first_gapped: uppercase seq; if first_gapped: profile
        alignment[0].push_back(aligned_seq[0]);
        // if not first_gapped: lowercase seq; if first_gapped: sequence
        alignment[1].push_back(aligned_seq[1]);
      }
    }
    //update number of performed alignments
    prev_alignments = next_alignments;
  }
  else {
    alignment = previous_alignment;
  }
  return alignment;
}


int msa::count_alignments(double identity_cutoff,
                          const std::vector<double>& identities) {
  auto count = 0;
  for (auto& item: identities) {
    if (item >= identity_cutoff) {
      ++count;
    }
  }
  return count;
}


std::vector<fasta::SequenceList> msa::perform_msa_round_gapped(
    const fasta::FastaData& fasta_data,
    const fasta::FastaData& fasta_data_alignment,
    const profile::ProfileMap& profile,
    const FeatureScores& f_profile,
    double gap_open_pen,
    double end_pen,
    double gap_ext_pen,
    double identity_cutoff,
    int codon_length,
    const std::vector<double>& identities,
    int& prev_alignments,
    const f_config::FeatureSettingsMap& f_set,
    std::vector<fasta::SequenceList> previous_alignment, int refine_seq,
    const bool no_feat)
{
  std::vector<fasta::SequenceList> alignment = {{}, {}};
  auto next_alignments = count_alignments(identity_cutoff, identities);
  if (next_alignments > prev_alignments) {
    fasta::SequenceList aligned_seq;
    for (size_t i = 0; i < fasta_data.sequences.size(); ++i) {
      if (identities[i] >= identity_cutoff) {
        // NW alignment of the ith seq against the profile
        aligned_seq = msa::align_pairwise(fasta_data.sequences[i],
                                          profile, f_profile, gap_open_pen,
                                          end_pen, gap_ext_pen, codon_length,
                                          true, no_feat);
        // profile
        alignment[0].push_back(aligned_seq[0]);
        // sequence
        alignment[1].push_back(aligned_seq[1]);
      }
    }
    //update number of performed alignments
    prev_alignments = next_alignments;
    alignment = merge_alignments(alignment);
  }
  else {
    alignment = previous_alignment;
  }
  return alignment;
}


std::vector<fasta::SequenceList> msa::merge_alignments(
    const std::vector<fasta::SequenceList>& pairwise_alignments)
{
  std::vector<fasta::SequenceList> multi_alignment;
  /// TODO: change it so that changing the format is no longer needed at the
  /// end of the function
  multi_alignment = {
          {pairwise_alignments[0][0], pairwise_alignments[0][0]},
          {pairwise_alignments[1][0], pairwise_alignments[1][0]}
  };

  assert(pairwise_alignments[0].size() == pairwise_alignments[1].size());
  for (size_t i = 1; i < pairwise_alignments[0].size(); ++i) {
    fasta::SequenceList pair_al = {pairwise_alignments[0][i],
                                   pairwise_alignments[1][i]};
    multi_alignment = add_alignment(multi_alignment, pair_al);
  }


  std::vector<fasta::SequenceList> result = {{}, {}};
  for (size_t i = 1; i < multi_alignment.size(); ++i) {
    result[0].push_back(multi_alignment[i][0]);
    result[1].push_back(multi_alignment[i][1]);
  }
  result = remove_gapcolumns(result);
  return result;
}


std::vector<fasta::SequenceList> msa::add_alignment(
    const std::vector<fasta::SequenceList>& multi_alignment,
    const fasta::SequenceList& pairwise_alignment)
{
  fasta::Sequence s;
  std::vector<fasta::SequenceList> merged(multi_alignment.size() + 1,
                                          fasta::SequenceList(2, s));
  auto i = 0;
  auto j = 0;
  const auto *profile1 = &multi_alignment[0][0].residues;
  const auto *profile2 = &pairwise_alignment[0].residues;
  int length1 = profile1->size();
  int length2 = profile2->size();

  int codon_length = profile1->at(0).codon.size();
  fasta::Residue gap_residue("-" + std::string(codon_length - 1, 'A'));
  while (i < length1 || j < length2)
  {
    if (i < length1 && j < length2
        && profile1->at(i).codon == profile2->at(j).codon)
    {
      for (size_t k = 0; k < multi_alignment.size(); ++k) {
        merged[k][0].residues.push_back(multi_alignment[k][0].residues[i]);
        merged[k][1].residues.push_back(multi_alignment[k][0].residues[i]);
      }
      merged[merged.size() - 1][0].residues.push_back(
          pairwise_alignment[1].residues[j]);
      merged[merged.size() - 1][1].residues.push_back(
          pairwise_alignment[1].residues[j]);
      ++i;
      ++j;
    } else if (i < length1 && profile1->at(i).codon[0] == '-') {
      merged[0][0].residues.push_back(gap_residue);
      merged[0][1].residues.push_back(gap_residue);
      merged[merged.size() - 1][0].residues.push_back(gap_residue);
      merged[merged.size() - 1][1].residues.push_back(gap_residue);
      for (size_t k = 1; k < multi_alignment.size(); ++k) {
        merged[k][0].residues.push_back(multi_alignment[k][0].residues[i]);
        merged[k][1].residues.push_back(multi_alignment[k][0].residues[i]);
      }
      ++i;
    } else if (j < length2 && profile2->at(j).codon[0] == '-') {
      merged[0][0].residues.push_back(gap_residue);
      merged[0][1].residues.push_back(gap_residue);
      merged[merged.size() - 1][0].residues.push_back(
          pairwise_alignment[1].residues[j]);
      merged[merged.size() - 1][1].residues.push_back(
          pairwise_alignment[1].residues[j]);
      for (size_t k = 1; k < multi_alignment.size(); ++k) {
        merged[k][0].residues.push_back(gap_residue);
        merged[k][1].residues.push_back(gap_residue);
      }
      ++j;
    }
  }
  return merged;
}


std::vector<fasta::SequenceList> msa::remove_gapcolumns(
    std::vector<fasta::SequenceList> alignment)
{
  std::vector<fasta::SequenceList> result = alignment;
  auto erased = 0;
  for (size_t i = 0; i < alignment[0][0].residues.size(); ++i) {
    if (alignment[0][0].residues[i].codon[0] == '-') {
      auto gaps = true;
      size_t j = 1;
      while (gaps && j < alignment[0].size()) {
        if (alignment[0][j].residues[i].codon[0] != '-') {
          gaps = false;
        }
        ++j;
      }
      if (gaps) {
        for (size_t k = 0; k < alignment[0].size(); ++k) {
          result[0][k].residues.erase(
              result[0][k].residues.begin() + i - erased);
          result[1][k].residues.erase(
              result[1][k].residues.begin() + i - erased);
        }
        ++erased;
      }
    }
  }
  return result;
}


bool msa::compare_alignments(const std::vector<fasta::SequenceList>& al1,
    const std::vector<fasta::SequenceList>& al2) {
  auto result = true;
  if (al1.size() != al2.size() || al1[0].size() != al2[0].size()) {
    result = false;
  }
  std::vector<std::string> al1_strings;
  std::vector<std::string> al2_strings;
  for (auto& item : al1) {
    for (auto& seq : item) {
      al1_strings.push_back(fasta::sequence_to_string(seq));
    }
  }
  for (auto& item : al2) {
    for (auto& seq : item) {
      al2_strings.push_back(fasta::sequence_to_string(seq));
    }
  }
  if (al1_strings != al2_strings) {
    result = false;
  }
  return result;
}
