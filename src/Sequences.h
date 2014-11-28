#ifndef SEQUENCES_H
#define SEQUENCES_H


#include "types.h"
#include <iostream>
#include <vector>


class Residue;
class Profile;
class FeaturesProfile;
class Sequences{
public:
	//constructor
  Sequences(codonSeqWithNamesList& s);
  Sequences();
	//getters
  seqNames get_names();
	//main functionality
	string_sequences performMSAfirstround(Profile& outputProfile, 
                                        FeaturesProfile& outputFeaturesProfile, 
                                        double penalty, double endPenalty, 
                                        double extensionPenalty, 
                                        bool weightsModeOn, int codon_length, 
                                        identitiesList& identities);

	void performMSAnextRounds(string_sequences* prevAlignment,
                            Profile& outputProfile, 
                            FeaturesProfile& outputFeaturesProfile, 
                            double penalty, double endPenalty, 
                            double extensionPenalty, 
                            bool weightsModeOn, double identityCutoff, 
                            int codon_length, identitiesList& identities, 
                            int& prev_alignments);
  void add_usr_features(rulesTuplesList& feature_rules);
private:
	//functions
	void removeGaps(sequence& alignmentWithLowercase, 
                  sequence& alignmentWithoutLowercase, 
                  sequenceList& alignment);
	void alignPairwise(sequence& alNoLower, 
                     sequence& alWithLower, 
                     sequence& seq2,
                     Profile& prf, FeaturesProfile& featPrf,
                     double penalty, double endPenalty, double extensionPenalty,
                     int deb, int codon_length);
	double calcIdentity(const sequence& alignedSequence);
	int countAlignments(double identity_cutoff, identitiesList& identities);
  void add_feature_indexes(FeaturesProfile& fprf);
	//variables 
	int m_seqNr;
	int m_firstSequenceSize;
	sequenceList m_sequences_aa;
	seqNames m_sequence_names;
};

#endif /* SEQUENCES_H */
