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
	std::vector<std::vector<double> > getMatrix();
	void printProfile();
	void createProfile(const std::vector< std::vector<std::string> >&,const std::vector<bool>&, const std::vector<double>&, bool);
	void expandListOfFeatures(const std::vector<std::string>&);
private:
	std::vector< std::vector<double> > prfMatrix;
	void countOccurences(const std::vector< std::vector<std::string> >&,const std::vector<bool>&);
	void countOccurences(const std::vector< std::vector<std::string> >&,const std::vector<double>&);
	int findFeaturesIndex(std::string);
	std::string name(std::string,int);
};
