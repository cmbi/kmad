#include <iostream>
#include <string>
#include <vector>
class Profile;
class ScoringMatrix{
public:
	ScoringMatrix(int,int,int); //constructor
	~ScoringMatrix();
	ScoringMatrix(ScoringMatrix&);
	ScoringMatrix operator=(ScoringMatrix&);
	void fillMatrix(std::string, std::string);
	void calculateScores(std::string,Profile&,int);
//	std::vector<std::string> 
	void nwAlignment(std::vector<std::string>*,std::string,Profile&, std::string);
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

