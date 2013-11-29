#include <iostream>
#include <string>
#include <vector>
class FeaturesProfile{
public:
	FeaturesProfile(std::vector< std::vector<double> >);
	FeaturesProfile();
	double getElement(int, std::string);
	double getElement(int, int);	
	double getScore(int,std::string);
	void printProfile();
	void createProfile(const std::vector< std::vector<std::string> >&,const std::vector<bool>&);
private:
	std::vector< std::vector<double> > prfMatrix;
	void countOccurences(const std::vector< std::vector<std::string> >&,const std::vector<bool>&);
	int findFeaturesIndex(std::string);
	void expandListOfFeatures(const std::vector< std::vector<std::string> >& , const std::vector<bool>&);
	std::string name(std::string,int);
};
