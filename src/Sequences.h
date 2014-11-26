#ifndef SEQUENCES_H
#define SEQUENCES_H
#include <iostream>
#include <vector>
class Residue;
class Profile;
class FeaturesProfile;
class Sequences{
public:
	//constructor
	Sequences(std::vector< std::vector< std::vector< std::string> > >& );
  Sequences();
	//getters
  std::vector<std::string> get_names();
	//main functionality
	std::vector<std::string> performMSAfirstround(Profile&, FeaturesProfile&, 
                                                double, double, double, 
                                                std::string, bool, int, 
                                                std::vector<double>&);
	void performMSAnextRounds(std::vector<std::string>*,
                            Profile&, FeaturesProfile&, 
                            double, double, double, 
                            std::string, bool, double, 
                            int, std::vector<double>&, 
                            int&);
	void printSequence(int) const;
	void add_usr_features(std::vector<std::tuple<std::string, std::string, int, 
                                               int, int, double, double, 
                                               double, double, std::string, 
                                               std::string > >&);
private:
	//functions
	void removeGaps(std::vector<Residue>&, std::vector<Residue>&, 
                  std::vector<std::vector<Residue>>&);
	void alignPairwise(std::vector<Residue>&, std::vector<Residue>&, 
                     std::vector<Residue>&,
                     Profile&, FeaturesProfile&,
                     double, double, double, int, std::string, int);
	double calcIdentity(const std::vector<Residue>&);
	int countAlignments(double, std::vector<double>&);
  void add_feature_indexes(FeaturesProfile&);
	//variables 
	int seqNr;
	int firstSequenceSize;
	std::vector<std::vector<Residue>> sequences_aa;
	std::vector<std::string> m_sequence_names;
};
#endif /* SEQUENCES_H */
