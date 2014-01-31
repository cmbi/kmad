#include <iostream>
#include <string>
#include <vector>
class Profile;
class FeaturesProfile;
class ScoringMatrix{
public:
	ScoringMatrix(int,int,double); //constructor
	ScoringMatrix(int,int,double,double); //constructor
	~ScoringMatrix();
	ScoringMatrix(ScoringMatrix&);
	ScoringMatrix operator=(ScoringMatrix&);
	void fillMatrix(std::string, std::string);
	void calculateScores(std::string,Profile&,int);
	void calculateScores(std::vector<std::string>,Profile&,FeaturesProfile&,int);
	void nwAlignment(std::vector<std::string>*,std::string,Profile&, std::string);
	void nwAlignment(std::vector<std::vector<std::string> >*,std::vector<std::string>,Profile&,FeaturesProfile&,std::string);
	std::vector<int> findBestScore();
	std::vector< std::vector<double> > getVec(); //getter
private:
	int iLength;
	int jLength;
	double gapOpening;
	double gapExtension;
	double gapOpeningHorizontal;
	double gapExtensionHorizontal;
	std::vector< std::vector<double> > matrixV,matrixG,matrixH;
};

