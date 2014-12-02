#ifndef SCORINGMATRIX_H
#define SCORINGMATRIX_H

#include "types.h"
#include <iostream>
#include <string>
#include <vector>

class Residue;
class Profile;
class FeaturesProfile;

typedef std::vector<double> ScoringMatrixRow;
typedef std::vector<ScoringMatrixRow> SingleScoringMatrix;
typedef std::vector<int> ValueCoords;

class ScoringMatrix{
public:
  /// 
  /// Constructor; creates scoring matrices for the two sequences
  /// (sequence and a profile)
  /// @param s1size length of the pseudo sequence (=profile)
  /// @param s2size lengtn of the sequence that will be aligned to the profile
  /// @param pen gap opening penalty
  /// @param extensionPenalty gap extension penalty
  /// @param endPenalty penalty for gaps at the beginning and the end
  ///
  ScoringMatrix(int s1_size, int s2_size, double gap_open_pen, double end_pen, 
                double gap_ext_pen);
  ///
  /// Fills in the scoring matrices m_matrix_v, m_matrix_g, m_matrix_h
  ///
  void CalculateScores(ResidueSequence s2, Profile& prf, 
                       FeaturesProfile& feat_prf, int codon_length);
  ///
  /// traces back the alignment path in the scoring matrices
  ///
  void PerformNWAlignment(SequenceList& result, ResidueSequence s2, 
                          Profile& prf, FeaturesProfile& feat_prf, 
                          int codon_length);
private:
  ///
  /// finds the best score either in the last column or in the last row of the 
  /// V matrix (takes the end gap penaltie into account)
  ///
  ValueCoords FindBestScore();
  int m_i_length;
  int m_j_length;
  double m_gap_opening;
  double m_gap_extension;
  double m_end_pen;
  double m_gap_opening_horizontal;
  double m_gap_extension_horizontal;
  SingleScoringMatrix m_matrix_v,m_matrix_g,m_matrix_h;
};

#endif /* SCORINGMATRIX_H */
