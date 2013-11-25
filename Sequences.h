#include <iostream>
#include <vector>
class Profile;
class Sequences{
public:
	//constructor
	Sequences(std::vector< std::vector<std::string> >);
	Sequences(std::vector< std::vector< std::vector<std::string> > >);
	Sequences();
	//
	void combinePairwiseAlignments(Sequences);  //MODIFY
	//getters
	std::vector< std::vector<std::string> > getSequences();
	std::vector<std::string> performMSA(Profile*,int,bool);
	void printEncodedSequence(int);
private:
	//functions
	void removeGaps(std::string *,std::string *,std::vector<std::string> &);
	void alignPairwise(std::string *,std::string *, std::string, Profile&, int, int,bool);
	std::vector< std::vector<std::string> > sequences;	
	std::vector< std::vector< std::vector<std::string> > > sequencesEncoded;	
	double calcIdentity(const std::string&);
	double countIdenticalResidues(std::vector<std::string>&);
	//variables 
	int seqNr;
	int firstSequenceSize;
};
