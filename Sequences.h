#include <iostream>
#include <vector>
class Profile;
class FeaturesProfile;
class Sequences{
public:
	//constructor
	Sequences(std::vector< std::vector< std::vector<std::string> > >);
	//getters
	std::vector< std::vector<std::string> > getSequences();
	std::vector< std::vector< std::vector<std::string> > > getEncodedSequences();
	//main functionality
	std::vector<std::string> performMSAencoded(Profile&,FeaturesProfile&,double,double,double,std::string, bool,int, std::vector<double>&);
	std::vector<std::string> performMSAnextRound(Profile&,FeaturesProfile&,double,double,double,std::string, bool, double,int, std::vector<double>);

private:
	//functions
	void removeGaps(std::vector<std::string> &,std::vector<std::string> &,std::vector<std::vector<std::string> >&);
	void alignPairwise(std::vector<std::string> &,std::vector<std::string> &, std::vector<std::string>, Profile&, FeaturesProfile&,double,double,double, int,std::string,int);
	std::vector< std::vector<std::string> > sequences;	
	std::vector< std::vector< std::vector<std::string> > > sequencesEncoded;	
	double calcIdentity(const std::vector<std::string>&);
	double countIdenticalResidues(std::vector<std::string>&);
	//variables 
	int seqNr;
	int firstSequenceSize;
	
};
