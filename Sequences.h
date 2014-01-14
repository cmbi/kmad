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
	~Sequences();
	Sequences operator=(Sequences&);
	//
	void combinePairwiseAlignments(Sequences);  //MODIFY
	//getters
	std::vector< std::vector<std::string> > getSequences();
	std::vector< std::vector< std::vector<std::string> > > getEncodedSequences();
	std::vector<std::string> performMSA(Profile*,int,std::string);
	std::vector<std::string> performMSAencoded(std::vector<std::vector<double> >*,std::vector<std::vector<double> >*,int,double,std::string, bool);
	std::vector<std::string> performMSAnextRound(Profile*,FeaturesProfile*,int,double,std::string, bool, double);
private:
	//functions
	void removeGaps(std::string *,std::string *,std::vector<std::string> &);
	void removeGaps(std::vector<std::string> *,std::vector<std::string> *,std::vector<std::vector<std::string> >&);
	void alignPairwise(std::string *,std::string *, std::string, Profile&, int, int,std::string);
	void alignPairwise(std::vector<std::string> *,std::vector<std::string> *, std::vector<std::string>, Profile&, FeaturesProfile&, int,double, int,std::string);
	std::vector< std::vector<std::string> > sequences;	
	std::vector< std::vector< std::vector<std::string> > > sequencesEncoded;	
	double calcIdentity(const std::string&);
	double calcIdentity(const std::vector<std::string>&);
	double countIdenticalResidues(std::vector<std::string>&);
	//variables 
	int seqNr;
	int firstSequenceSize;
	
};
