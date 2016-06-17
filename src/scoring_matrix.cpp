#include "compare_doubles.h"
#include "feature_scores.h"
#include "profile.h"
#include "scoring_matrix.h"

#include <boost/filesystem.hpp>


ScoringMatrix::ScoringMatrix(int s1_size,int s2_size, double pen,
                             double end_pen, double gap_ext_pen,
                             const bool no_feat)
: m_i_length(s1_size),
  m_j_length(s2_size),
  m_gap_opening(pen),
  m_gap_extension(gap_ext_pen),
  m_end_pen(end_pen),
  m_no_feat(no_feat)
{
  //creates a row for the scoring matrices of length m_j_length
  //(length of the jth sequence + 1)
  ScoringMatrixRow row(m_j_length + 1,0);
  //creates a vector of vectors 'row', of length m_i_length+1
  //(length of the ith sequence +1)
  m_matrix_v.assign(m_i_length + 1, row);
  m_matrix_g.assign(m_i_length + 1, row);
  m_matrix_h.assign(m_i_length + 1, row);
}


void ScoringMatrix::calculate_scores(const fasta::Sequence& sequence,
                                     const profile::ProfileMap& profile,
                                     const FeatureScores& f_profile,
                                     int codon_length) {
  assert(m_matrix_v.size() == m_matrix_g.size());
  assert(m_matrix_v.size() == m_matrix_h.size());

  for (unsigned int i = 1; i < m_matrix_v.size(); ++i) {
    m_matrix_v[i][0] = -10000000; //=== -infinity
    m_matrix_h[i][0] = i * m_end_pen;
    m_matrix_g[i][0] = -10000000;
  }
  for (unsigned int i = 1; i < m_matrix_v[0].size(); ++i) {
    m_matrix_v[0][i] = -10000000;
    m_matrix_h[0][i] = -10000000;
    m_matrix_g[0][i] = i * m_end_pen;
  }

  double score1, score2, score3 = 0;
  for (size_t i = 1; i < m_matrix_v.size(); ++i) {
    for (size_t j = 1; j < m_matrix_v[i].size(); ++j) {
      ///V
      // double profile_score =
      //   profile.at(sequence.residues[j - 1].codon[0])[i - 1];
      double profile_score = profile::get_score(profile, i - 1,
          sequence.residues[j - 1].codon[0]);
      double feature_score = 0;
      if (!m_no_feat) {
        for (auto& feat_name : sequence.residues[j - 1].features) {
          feature_score += f_profile.get_score(feat_name, i - 1);
        }
      }

      score1 = m_matrix_v[i-1][j-1];
      score2 = m_matrix_g[i-1][j-1];
      score3 = m_matrix_h[i-1][j-1];

      m_matrix_v[i][j] = std::max<double>(score1,
        std::max<double>(score2, score3)) + profile_score + feature_score;
      ///G
      score1 = m_matrix_v[i-1][j] + m_gap_opening;
      score2 = m_matrix_h[i-1][j] + m_gap_extension;
      score3 = m_matrix_g[i-1][j] + m_gap_opening;
      m_matrix_h[i][j] = std::max<double>(score1, std::max<double>(score2,
                                                                   score3));
      ///H
      score1 = m_matrix_v[i][j-1] + m_gap_opening;
      score2 = m_matrix_g[i][j-1] + m_gap_extension;
      score3 = m_matrix_h[i][j-1] + m_gap_opening;
      m_matrix_g[i][j] = std::max<double>(score1, std::max<double>(score2,
                                                                   score3));
    }
  }
}


ValueCoords ScoringMatrix::find_best_score() {
  int max_i = m_matrix_v.size()-1;
  int max_j = m_matrix_v[0].size()-1;
  int n = max_i; // last row of m_matrix_v
  int m = max_j; // last column of m_matrix_v
  double max_i_val = m_matrix_v[max_i][max_j];
  double max_j_val = m_matrix_v[max_i][max_j];
  double real_val;
  double max = m_matrix_v[max_i][max_j];
  for (int i = 0; i <= n ; ++i) {    //finds max score in the last row
    // add end gap penalties to the score to calc the 'real' score
    // of the alignment
    real_val = m_matrix_v[i][m] + m_end_pen * (m_matrix_v.size() - i - 1);
    if (real_val > max) {
      max_i_val = real_val;
      max = real_val;
      max_i = i;
    }
  }
  for (int i = 0; i <= m; ++i) {    //finds max score in the last column
    real_val = m_matrix_v[n][i] + m_end_pen * (m_matrix_v[0].size() - i - 1);
    // add end gap penalties to the score, to calc the 'real' score
    // of the alignment
    if (real_val > max) {
      max_j_val = real_val;
      max = real_val;
      max_j = i;
    }
  }
  ValueCoords res_arr;
  if (max_i_val > max_j_val) {      //max score is in the last row
    res_arr.push_back(max_i);
    res_arr.push_back(m);
  } else {          //max score is in the last column
    res_arr.push_back(n);    // max score position (i)
    res_arr.push_back(max_j);    //max score position (j)
  }
  return res_arr; //coords of the max score
}

fasta::SequenceList ScoringMatrix::backtrace_alignment_path(
    const fasta::Sequence& sequence, const profile::ProfileMap& profile,
    const FeatureScores& f_profile,
    int codon_length) {
  //creating polyA pseudoSequence representing the profile,
  //to know later where are the gaps in the profile
  fasta::Residue ala('A' + std::string(codon_length - 1, 'A'), {});
  fasta::Residue gap_residue('-' + std::string(codon_length - 1, 'A'), {});

  fasta::Sequence new_s1;
  fasta::Sequence new_s2;
  const size_t profile_length = profile.begin()->second.size();
  // TODO: i and j aren't great names...are they?
  size_t i = profile_length;
  size_t j = sequence.residues.size();

  //if bestScore isn't in the lower right corner, then add gaps
  //to new_s1 or new_s2
  ValueCoords best_score = find_best_score();

  // TODO: comparing value to size of something? fishy!
  if (best_score[0] != (signed)m_matrix_v.size()-1
      || best_score[1] != (signed)m_matrix_v[0].size()-1) {
    i = best_score[0];
    j = best_score[1];
    for (size_t k = profile_length; k > i; --k) {
      new_s1.residues.push_back(ala);
      new_s2.residues.push_back(gap_residue);
    }
    for (size_t k = sequence.residues.size(); k > j; --k) {
      new_s2.residues.push_back(sequence.residues[k - 1]);
      new_s1.residues.push_back(gap_residue);
    }
  }
  std::string current_matrix = "V";

  assert(m_matrix_v.size() == m_matrix_g.size());
  assert(m_matrix_v.size() == m_matrix_h.size());
  fasta::Residue new_res1;
  fasta::Residue new_res2;

  //trace back the matrix
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && current_matrix == "V") {  //match/mismatch
      new_res1 = ala;
      new_res2 = sequence.residues[j - 1];
      // double profile_score = profile.at(
      //     sequence.residues[j - 1].codon[0])[i - 1];
      double profile_score = profile::get_score(profile, i - 1,
          sequence.residues[j - 1].codon[0]);
      double feature_score = 0;
      if (!m_no_feat) {
        for (auto& feat_name : sequence.residues[j - 1].features) {
          feature_score += f_profile.get_score(feat_name, i - 1);
        }
      }
      double final_score = profile_score + feature_score;
      if (m_matrix_v[i][j] != m_matrix_v[i-1][j-1] + final_score) {
        if (i > 0 && j > 0
            && compare_doubles::is_equal(m_matrix_v[i][j],
                                         m_matrix_g[i-1][j-1] + final_score)) {
          current_matrix = "G";
        } else if (i > 0 && j > 0
                   && compare_doubles::is_equal(m_matrix_v[i][j],
                                                m_matrix_h[i-1][j-1]
                                                + final_score)) {
          current_matrix = "H";
        }
      }
      --i;
      --j;
    } else if (j > 0 && current_matrix == "G") {  //gap in seq2
      new_res1 = gap_residue;
      new_res2 = sequence.residues[j - 1];
      if (compare_doubles::is_equal(m_matrix_g[i][j],
                                    m_matrix_v[i][j - 1] + m_gap_opening)) {
        current_matrix = "V";
      }
      else if (compare_doubles::is_equal(m_matrix_g[i][j],
                                         m_matrix_h[i][j - 1]
                                         + m_gap_opening)) {
        current_matrix = "H";
      }
      --j;
    } else if (i > 0 && current_matrix == "H") {  //gap in profile
      new_res1 = ala;
      new_res2 = gap_residue;
      if (compare_doubles::is_equal(m_matrix_h[i][j],
                                    m_matrix_v[i - 1][j] + m_gap_opening)) {
        current_matrix = "V";
      }
      else if (compare_doubles::is_equal(m_matrix_h[i][j],
                                         m_matrix_g[i - 1][j]
                                         + m_gap_opening)) {
        current_matrix = "G";
      }
      --i;
    }
    new_s1.residues.push_back(new_res1);
    new_s2.residues.push_back(new_res2);
  }
  //need to reverse the sequences, because tracing back the alignment goes
  //from the end to the beginning
  std::reverse(new_s1.residues.begin(), new_s1.residues.end());
  std::reverse(new_s2.residues.begin(), new_s2.residues.end());
  return {new_s1, new_s2};
}

std::vector<ScoringMatrixRow> ScoringMatrix::get_V_matrix() {
  return m_matrix_v;
}
