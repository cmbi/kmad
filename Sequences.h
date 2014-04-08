#include <iostream>
#include <vector>
class Profile;
class FeaturesProfile;
class Sequences{
public:
	//constructor, copy constructor, destructor, operator=
	Sequences(std::vector< std::vector<std::string> >);
	Sequences(std::vector< std::vector< std::vector<std::string> > >);
	Sequences(Sequences&);
	Sequences();
	//getters
	std::vector< std::vector<std::string> > getSequences();
	std::vector< std::vector< std::vector<std::string> > > getEncodedSequences();
	//main functionality
	std::vector<std::string> performMSAencoded(Profile&,FeaturesProfile&,double,double,std::string, bool,int,int,int, std::vector<double>&);
	std::vector<std::string> performMSAnextRound(Profile&,FeaturesProfile&,double,double,std::string, bool, double,int,int,int, std::vector<double>);

private:
	//functions
	void removeGaps(std::vector<std::string> &,std::vector<std::string> &,std::vector<std::vector<std::string> >&);
	void alignPairwise(std::vector<std::string> &,std::vector<std::string> &, std::vector<std::string>, Profile&, FeaturesProfile&, double,double, int,std::string,int);
	std::vector< std::vector<std::string> > sequences;	
	std::vector< std::vector< std::vector<std::string> > > sequencesEncoded;	
	double calcIdentity(const std::vector<std::string>&);
	double countIdenticalResidues(std::vector<std::string>&);
	//variables 
	int seqNr;
	int firstSequenceSize;
	
};
