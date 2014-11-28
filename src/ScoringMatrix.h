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
  // Constructor
  //
  // @param int int1
	ScoringMatrix(int s1size, int s2size, double pen, double endPenalty, 
                double extensionPenalty); //constructor
	void calculateScores(sequence s2, Profile& prf, FeaturesProfile& featPrf, 
                       int debug, int codon_length);
	void nwAlignment(sequenceList* result, sequence s2, 
                   Profile& prf, FeaturesProfile& featPrf, int codon_length);
private:
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
