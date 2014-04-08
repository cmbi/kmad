#include <iostream>
#include <string>
#include <vector>
class Profile;
class FeaturesProfile;
class ScoringMatrix{
public:
	ScoringMatrix(int,int,double,double); //constructor
	void calculateScores(std::vector<std::string>,Profile&,FeaturesProfile&,int, int);
	void nwAlignment(std::vector<std::vector<std::string> >*,std::vector<std::string>,Profile&,FeaturesProfile&,std::string,int);
	std::vector< std::vector<double> > getVec(); //getter
private:
	int iLength;
	int jLength;
	double gapOpening;
	double gapExtension;
	double gapOpeningHorizontal;
	double gapExtensionHorizontal;
	std::vector< std::vector<double> > matrixV,matrixG,matrixH;
	std::vector<int> findBestScore();
};

