#include "features_profile.h"
#include "profile.h"
#include "scoring_matrix.h"
#include "substitution_matrix.h"
#include "txtproc.h"
#include "vec_util.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


ScoringMatrix::ScoringMatrix(int s1_size,int s2_size, double pen,
                             double end_pen, double gap_ext_pen)
: m_i_length(s1_size),
  m_j_length(s2_size),
  m_gap_opening(pen),
  m_gap_extension(gap_ext_pen),
  m_end_pen(end_pen)
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


// TODO: Implement
void ScoringMatrix::calculate_scores(const fasta::Sequence& sequence,
                                     const ProfileMap& profile,
                                     const FeaturesProfile& f_profile,
                                     int codon_length) {
  assert(m_matrix_v.size() == m_matrix_g.size());
  assert(m_matrix_v.size() == m_matrix_h.size());

  for (unsigned int i = 1; i < m_matrix_v.size(); i++) {
    m_matrix_v[i][0] = -10000000; //=== infinity
    m_matrix_h[i][0] = -10000000;
    m_matrix_g[i][0] = i * m_end_pen;
  }
  for (unsigned int i = 1; i < m_matrix_v[0].size(); i++) {
    m_matrix_v[0][i] = -10000000;
    m_matrix_h[0][i] = i * m_end_pen;
    m_matrix_g[0][i] = -10000000;
  }
  double score1, score2, score3 = 0;
  for (size_t i = 1; i < m_matrix_v.size(); i++) {
    for (size_t j = 1; j < m_matrix_v[i].size(); j++) {
      ///V
      //double prf_score = prf.get_element(i - 1,
      //                                   sequence.residues[j - 1].codon[0]);
      double profile_score = 
        profile.at(sequence.residues[j-1].codon[0])[i - 1];
      double feature_score = 0;
      for (auto& feat_name : sequence.residues[j - 1].features) {
        feature_score += f_profile.m_scores.at(feat_name)[i - 1]; 
      }
      // double add_score = 0;
      // FeaturesList features = s2[j].get_feat_indexes();
      // feat_prf.get_score(i-1, features, add_score);

      double final_score = profile_score + feature_score;
      score1 = m_matrix_v[i-1][j-1];
      score2 = m_matrix_g[i-1][j-1];
      score3 = m_matrix_h[i-1][j-1];

      m_matrix_v[i][j] = std::max<double>(score1,
        std::max<double>(score2, score3)) + final_score;
      ///G
      score1 = m_matrix_v[i-1][j] + m_gap_opening;
      score2 = m_matrix_g[i-1][j] + m_gap_extension;
      m_matrix_g[i][j] = (score1 > score2) ? score1 : score2;
      ///H
      score1 = m_matrix_v[i][j-1] + m_gap_opening;
      score2 = m_matrix_h[i][j-1] + m_gap_extension;
      m_matrix_h[i][j] = (score1 > score2) ? score1 : score2;
    }
  }
}


// ValueCoords ScoringMatrix::FindBestScore() {
//   int max_i = m_matrix_v.size()-1;
//   int max_j = m_matrix_v[0].size()-1;
//   int n = max_i; // last row of m_matrix_v
//   int m = max_j; // last column of m_matrix_v
//   double max_i_val = m_matrix_v[max_i][max_j];
//   double max_j_val = m_matrix_v[max_i][max_j];
//   double real_val;
//   double max = m_matrix_v[max_i][max_j];
//   for (int i = 0; i < n ; i++) {    //finds max score in the last row
//     // add end gap penalties to the score to calc the 'real' score
//     // of the alignment
//     real_val = m_matrix_v[i][m] + m_end_pen * (m_matrix_v.size() - i);
//     if (real_val > max) {
//       max_i_val = real_val;
//       max = real_val;
//       max_i = i;
//     }
//   }
//   for (int i = 0; i < m; i++) {    //finds max score in the last column
//     real_val = m_matrix_v[n][i] + m_end_pen * (m_matrix_v[0].size() - i);
//     // add end gap penalties to the score, to calc the 'real' score
//     // of the alignment
//     if (real_val > max) {
//       max_j_val = real_val;
//       max = real_val;
//       max_j = i;
//     }
//   }
//   ValueCoords res_arr;
//   if (max_i_val > max_j_val) {      //max score is in the last row
//     res_arr.push_back(max_i);
//     res_arr.push_back(m);
//   } else {          //max score is in the last column
//     res_arr.push_back(n);    // max score position (i)
//     res_arr.push_back(max_j);    //max score position (j)
//   }
//   return res_arr; //coords of the max score
// }
// 
// // TODO: Uncomment. This is required.
// void ScoringMatrix::PerformNWAlignment(std::vector<fasta::Sequence>& result,
//                                        fasta::Sequence s2, Profile& prf,
//                                        FeaturesProfile& feat_prf,
//                                       int codon_length) {
// /*
//   //creating polyA pseudoSequence representing the profile,
//   //to know later where are the gaps in the profile
//   Sequence s1(prf.get_matrix()[0].size()+1, Residue('A', codon_length));
//   Residue gap_residue = Residue('-', codon_length);
//   fasta::Residue gap_residue('-' + std::string(codon_length, 'A'), {});
//   s2.insert(s2.begin(), gap_residue);
//   Sequence new_s1;
//   Sequence new_s2;
//   SequenceList ali; //alignment
//   Residue new_char1;
//   Residue new_char2;
//   int i = s1.size()-1;
//   int j = s2.size()-1;
//   std::string current_matrix = "V";
//   //if bestScore isn't in the lower right corner, then add gaps
//   //to new_s1 or new_s2
//   ValueCoords best_score = FindBestScore();
//   if (best_score[0] != (signed)m_matrix_v.size()-1
//       || best_score[1] != (signed)m_matrix_v[0].size()-1) {
//     i = best_score[0];
//     j = best_score[1];
//     for (int k = s1.size()-1; k > i; k--) {
//       new_char1 = s1[k];
//       new_char2 = gap_residue;
//       new_s1.push_back(new_char1);
//       new_s2.push_back(new_char2);
//     }
//     for (int k = s2.size()-1; k > j; k--) {
//       new_char1 = gap_residue;
//       new_char2 = s2[k];
//       new_s2.push_back(new_char2);
//       new_s1.push_back(new_char1);
//     }
//   }
// 
//   assert(m_matrix_v.size() == m_matrix_g.size());
//   assert(m_matrix_v.size() == m_matrix_h.size());
// 
//   //trace back the matrix
//   while (i > 0 || j > 0) {
//     if (i > 0 && j > 0 && current_matrix == "V") {  //match/mismatch
//       new_char1 = s1[i];
//       new_char2 = s2[j];
//       double prf_score = prf.get_element(i-1,s2[j].get_aa());
//       double add_score = 0;
//       FeaturesList features = s2[j].get_feat_indexes();
//       feat_prf.get_score(i-1, features, add_score);
//       double final_score = prf_score + add_score;
//       if (m_matrix_v[i][j] != m_matrix_v[i-1][j-1] + final_score) {
//         if (i > 0 && j > 0
//             && m_matrix_v[i][j] == m_matrix_g[i-1][j-1]+final_score) {
//           current_matrix = "G";
//         } else if (i > 0 && j > 0
//                    && m_matrix_v[i][j] == m_matrix_h[i-1][j-1]+final_score) {
//           current_matrix = "H";
//         }
//       }
//       i--;
//       j--;
//     } else if (i > 0 && current_matrix == "G") {  //gap in seq2
//       new_char1 = s1[i];
//       new_char2 = gap_residue;
//       if (m_matrix_g[i][j] == m_matrix_v[i-1][j] + m_gap_opening)
//         current_matrix = "V";
//       i--;
//     } else if (j > 0 && current_matrix == "H") {  //gap in profile
//       new_char1 = gap_residue;
//       new_char2 = s2[j];
//       if (m_matrix_h[i][j] == m_matrix_v[i][j-1] + m_gap_opening_horizontal) {
//         current_matrix = "V";
//       }
//       j--;
//     }
//     new_s1.push_back(new_char1);
//     new_s2.push_back(new_char2);
//   }
//   std::reverse(new_s1.begin(),new_s1.end()); //need to reverse the sequences, because tracing back the alignment goes from the end to the beginning
//   std::reverse(new_s2.begin(),new_s2.end());
//   ali.push_back(new_s1);
//   ali.push_back(new_s2);
//   result = ali;
//   */
// }
