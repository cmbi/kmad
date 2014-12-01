#ifndef SCORINGMATRIX_H
#define SCORINGMATRIX_H

#include "types.h"
#include <iostream>
#include <string>
#include <vector>
class Residue;
class Profile;
class FeaturesProfile;
typedef std::vector<double> scoringMatrixRow;
typedef std::vector<scoringMatrixRow> scoringMatrix;
typedef std::vector<int> valueCoords;
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
	ScoringMatrix(int s1size, int s2size, double pen, double endPenalty, 
                double extensionPenalty);
  ///
  /// Fills in the scoring matrices m_matrixV, m_matrixG, m_matrixH
  ///
	void calculateScores(ResidueSequence s2, Profile& prf, 
                       FeaturesProfile& featPrf, int codon_length);
  ///
  /// traces back the alignment path in the scoring matrices
  ///
	void nwAlignment(SequenceList* result, ResidueSequence s2, 
                   Profile& prf, FeaturesProfile& featPrf, int codon_length);
private:
  ///
  /// finds the best score either in the last column or in the last row of the 
  /// V matrix (takes the end gap penaltie into account)
  ///
	valueCoords findBestScore();
	int m_iLength;
	int m_jLength;
	double m_gapOpening;
	double m_gapExtension;
	double m_endGapPenalty;
	double m_gapOpeningHorizontal;
	double m_gapExtensionHorizontal;
	scoringMatrix m_matrixV,m_matrixG,m_matrixH;
};

#endif /* SCORINGMATRIX_H */
