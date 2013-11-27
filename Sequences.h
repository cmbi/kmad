#include <iostream>
#include <vector>
class Profile;
class Sequences{
public:
	//constructor, copy constructor, destructor, operator=
	Sequences(std::vector< std::vector<std::string> >);
	Sequences(Sequences&);
	~Sequences();
	Sequences operator=(Sequences&);
	//
	void combinePairwiseAlignments(Sequences);  //MODIFY
	//getters
	std::vector< std::vector<std::string> > getSequences();
	std::vector<std::string> performMSA(Profile*,int,bool);
private:
	//functions
	void removeGaps(std::string *,std::string *,std::vector<std::string> &);
	void alignPairwise(std::string *,std::string *, std::string, Profile&, int, int,bool);
	std::vector< std::vector<std::string> > sequences;	
	double calcIdentity(const std::string&);
	double countIdenticalResidues(std::vector<std::string>&);
	//variables 
	int seqNr;
	int firstSequenceSize;
};
