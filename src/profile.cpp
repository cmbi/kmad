#include "profile.h"

#include "substitution_matrix.h"
#include "vec_util.h"

#include <boost/range/numeric.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


namespace {
  static const std::vector<char> ALPHABET = {'A','R','N','D','C','Q','E','G',
                                             'H','I','L','K','M','F','P','S',
                                             'T','W','Y','V'};
}

ProfileMap create_score_profile(const fasta::SequenceList& sequences) {
  ProfileMap p = create_profile(sequences);

  // Convert the profile occurrences to probabilities.
  for (auto& occ: p) {
    size_t i = 0;
    for (auto& v: occ.second) {
      occ.second[i] = v / sequences.size();
      i++;
    }
  }
  for (unsigned i = 0; i < p.size(); i++) {
    for (unsigned j = 0; j < p[i].size(); j++) {
      p[i][j] = 5;
    }
  }
  // for (unsigned i = 0; i < p[0].size(); i++) {
  //   std::vector<double> score_column(20, 0); 
  //   for (auto &prob: p){
  //     std::vector<double> sbst_column = substitution_matrix::get_column(
  //         prob.first);
  //     std::for_each(sbst_column.begin(), sbst_column.end(), [&](double &val){
  //           val *= prob.second[i];
  //         });
  //     for (size_t k = 0; k < sbst_column.size(); ++k) {
  //       score_column[k] += sbst_column[k] * prob.second[i];
  //     }
  //   }
  //   // replace the old column (with probs) with the new column (with scores)
  //   for (size_t j = 0; j < p.size(); ++j) {
  //     //p[j][i] = score_column[j];
  //     p[j][i] = 5;
  //   }
  // }
  // TODO: Calculate scores for each probability

  return p;
}


ProfileMap create_profile(const fasta::SequenceList& sequences) {
  ProfileMap p;

  // Initialise profile map with all valid letters to a vector of the correct
  // size, with all values set to 0.
  for (auto& c: ALPHABET) {
    p[c] = std::vector<double>(sequences[0].residues.size(), 0);
  }

  for (size_t i = 0; i < sequences[0].residues.size(); i++) {
    for (size_t j = 0; j < sequences.size(); j++) {
      char amino_acid = sequences[j].residues[i].codon[0];

      if (amino_acid == '-') {
        continue;
      }

      if (amino_acid == 'B') {
        assert(p.find('D') != p.end() && p.find('N') != p.end());
        p['D'][i] += 0.5;
        p['N'][i] += 0.5;
      } else if (amino_acid == 'Z') {
        assert(p.find('E') != p.end() && p.find('G') != p.end());
        p['E'][i] += 0.5;
        p['G'][i] += 0.5;
      } else if (amino_acid == 'X') {
        for (auto& kv: p) {
          kv.second[i] += 0.05;
        }
      } else {
        assert(p.find(amino_acid) != p.end());
        p[amino_acid][i] += 1.0;
      }
    }
  }

  return p;
}


Profile::Profile(ProfileMatrix mat)
:  m_prf_matrix(mat) {
}


Profile::Profile() {
}


ProfileMatrix Profile::get_matrix() const{
  return m_prf_matrix;
}


void Profile::ProcessProfile(std::vector<fasta::Sequence>& alignment) {
  CreateProfile(alignment);
  ProfileMatrix new_profile;
  for (unsigned int i = 0; i < m_prf_matrix[0].size(); i++) {
    Matrix2D columns_to_add;
    for (unsigned int j = 0; j < m_prf_matrix.size(); j++) {
      if (m_prf_matrix[j][i] != 0) {
        SbstMatColumn column_int;
        substitution_matrix::get_column(j, column_int);
        ProfileMatrixColumn column_j = vec_util::ConvertIntVecToDoubleVec(
            column_int);
        vec_util::MultiplyVectorByAScalar(column_j, m_prf_matrix[j][i]);
        columns_to_add.push_back(column_j);
      }
    }
    //add up columns from substitution matrix for amino acids seen on ith
    //position(times occurence/totalNrOfSeq))
    new_profile.push_back(vec_util::AddUp(columns_to_add));
  }
  vec_util::TransposeVec(new_profile);
  m_prf_matrix = new_profile;
}


// TODO: Do this in the constructor or use factory method
// TODO: Why is the parameter called 'alignment'?
// TODO: Implement
void Profile::CreateProfile(std::vector<fasta::Sequence>& alignment) {
  //ProfileMatrix tmp_result;
  //int no_of_sequences = alignment.size();
  //for (unsigned int i = 0; i < alignment[0].size(); i++) {
    //ProfileMatrixColumn profile_column(20,0);
    //int non_gaps = 0;
    //for (unsigned int j = 0; j < alignment.size(); j++) {
      //char seq_char(alignment[j][i].get_aa());
      //if (seq_char != '-') {
        ////either D or N, so add half a point to both
        //if (seq_char == 'B') {
          //profile_column[2] += 0.5;
          //profile_column[3] += 0.5;
        //} else if (seq_char == 'Z') {
          ////either D or N, so add half a point to both
          //profile_column[6] += 0.5;
          //profile_column[7] += 0.5;
        //} else if (seq_char == 'X') {
          //for (unsigned int k = 0; k < profile_column.size(); k++) {
            //profile_column[k] += 0.05;
          //}
        //} else {
          //int aacid_index = substitution_matrix::FindAminoAcidsIndex(seq_char);
          //profile_column[aacid_index] += 1;
        //}
        //non_gaps++;
      //}
    //}
    //vec_util::DivideVectorByAScalar(profile_column,no_of_sequences);
    ////vec_util::DivideVectorByAScalar(profileColumn,nonGaps);
    //tmp_result.push_back(profile_column);
  //}
  //m_prf_matrix = tmp_result;
  //vec_util::TransposeVec(m_prf_matrix);
}


double Profile::get_element(int position, char aacid) {
  double result;
  if (aacid == 'B') {
    // take half the score for asparagine and half the score for aspartate
    result = 0.5 * m_prf_matrix[2][position] + 0.5 * m_prf_matrix[3][position];
  } else if (aacid == 'Z') {
    // take half the score for glutamine and half the score for glutamate
    result = 0.5 * m_prf_matrix[6][position] + 0.5 * m_prf_matrix[7][position];
  } else if (aacid == 'X') {
    // take average score from scores for all residues
    result = 0;
    for (unsigned int i = 0; i < m_prf_matrix.size(); i++) {
      result += 0.05 * m_prf_matrix[i][position];
    }
  } else {
    // it's not any of the {B,Z,X} -> single amino acid
    int aacid_index = substitution_matrix::FindAminoAcidsIndex(aacid);
    result = m_prf_matrix[aacid_index][position];
  }
  return result;
}


double Profile::get_element(int aacid_index, int position) {
  return m_prf_matrix[aacid_index][position];
}
