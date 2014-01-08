#include <iostream>
#include <string>
#include <vector>
class Profile{
public:
	Profile(std::vector< std::vector<double> >); //constructor
	Profile();
	Profile(Profile&); //copy constructor
	Profile operator=(const Profile&);
//	Profile& operator=(Profile&);
	~Profile();
	void buildPseudoProfile(std::vector< std::string >&, const std::vector<bool>&);
	void buildPseudoProfile(std::vector< std::vector< std::string > >&, const std::vector<bool>&, const std::vector<double>&, bool);
	//getters//
	std::string getConsensusSequence();
	std::vector< std::vector<double> > getMatrix() const;
	double getElement(int, char);
	double getElement(int, int);
	void printProfile(int,int);
	void printProfile();
private:
	//functions
	void createProfile(std::vector<std::string>&,const std::vector<bool>&);
	void createProfile(std::vector< std::vector<std::string> >&,const std::vector<bool>&, const std::vector<double>&, bool);
//	std::vector< std::vector<double> > 
	void countOccurences(std::vector< std::vector<double> >&,std::vector<std::string>&,const std::vector<bool>&);
	void countOccurences(std::vector< std::vector<double> >&,std::vector< std::vector<std::string> >&,const std::vector<bool>&);
	void countOccurences(std::vector< std::vector<double> >&,std::vector< std::vector<std::string> >&,const std::vector<bool>&,const std::vector<double>&);
	double countNonGaps(int);
	int getMaxDoubleValue(std::vector<double>);
	//variables
	std::vector< std::vector<double> > prfMatrix;
};
