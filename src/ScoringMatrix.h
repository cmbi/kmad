#ifndef SCORINGMATRIX_H
#define SCORINGMATRIX_H

#include <iostream>
#include <string>
#include <vector>
class Residue;
class Profile;
class FeaturesProfile;
class ScoringMatrix{
public:
	ScoringMatrix(int, int, double, double, double); //constructor
	void calculateScores(std::vector<Residue>, Profile&, FeaturesProfile&, 
                       int, int);
	void nwAlignment(std::vector<std::vector<Residue>>*, std::vector<Residue>, 
                   Profile&, FeaturesProfile&, int);
	std::vector< std::vector<double> > getVec(); //getter
private:
	std::vector<int> findBestScore();
	int m_iLength;
	int m_jLength;
	double m_gapOpening;
	double m_gapExtension;
	double m_endGapPenalty;
	double m_gapOpeningHorizontal;
	double m_gapExtensionHorizontal;
	std::vector< std::vector<double> > m_matrixV,m_matrixG,m_matrixH;
};

#endif /* SCORINGMATRIX_H */
