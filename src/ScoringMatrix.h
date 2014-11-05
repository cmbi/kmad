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
                       int, int, int);
	void nwAlignment(std::vector<std::vector<Residue> >*, std::vector<Residue>, 
                   Profile&, FeaturesProfile&, std::string, int, int);
	std::vector< std::vector<double> > getVec(); //getter
private:
	int iLength;
	int jLength;
	double gapOpening;
	double gapExtension;
	double endGapPenalty;
	double gapOpeningHorizontal;
	double gapExtensionHorizontal;
	std::vector< std::vector<double> > matrixV,matrixG,matrixH;
	std::vector<int> findBestScore();
};

#endif /* SCORINGMATRIX_H */
