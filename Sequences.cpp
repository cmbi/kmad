#include "Sequences.h"
#include "FeaturesProfile.h"
#include "Profile.h"
#include "ScoringMatrix.h"
#include "substitutionMatrix.h"
#include "vecUtil.h"
#include <iostream>
#include <string>
#include <vector>
// 1st constructor - assigns s to vector< vector<string> > sequences - [[seqName,sequence],...]
Sequences::Sequences(std::vector< std::vector<std::string> > s){
	sequences = s;
	seqNr = s.size(); 
	firstSequenceSize = s.at(0).at(1).size();
}
//constructor for ENCODED SEQUENCES
Sequences::Sequences(std::vector< std::vector< std::vector<std::string> > > s){
	sequencesEncoded = s;
	seqNr = s.size(); 
	firstSequenceSize = s.at(0).at(1).size();
}
Sequences::Sequences(Sequences& that){
	sequences = that.sequences;
	sequencesEncoded = that.sequencesEncoded;
	seqNr = that.seqNr;
	firstSequenceSize = that.firstSequenceSize;
}
Sequences::Sequences(){}
Sequences::~Sequences(){
}
Sequences Sequences::operator=(Sequences& that){
	Sequences newSequences(that.sequences);
	newSequences.sequencesEncoded =  that.sequencesEncoded;
	newSequences.seqNr = that.seqNr;
	newSequences.firstSequenceSize = that.firstSequenceSize;
	return newSequences;
}
//function performMSA - returns ready alignment, with gaps cut out from the first sequence and so on...
std::vector<std::string> Sequences::performMSA(Profile& outputProfile,double penalty, double extensionPenalty, std::string verbose){
	outputProfile = Profile(substitutionMatrix::convertToProfileFormat(sequences.at(0).at(1))); //IMPLEMENT convertToProfileFormat // psuedoProfile - it's the one that will later be combined profile and sbst matrix
	std::vector<std::string> alignmentWithoutLowercase;	//working alignment - without lowercase around cut out residues - would make latter aligning more complicated
	std::vector<std::string> alignmentWithLowercase;		//lowercase before and after cut out residues -- final result 
	alignmentWithoutLowercase.push_back(sequences.at(0).at(1));
	alignmentWithLowercase.push_back(sequences.at(0).at(1));
	std::vector<bool> sequenceIdentity; 			//'true' stored for every sequence which identity with the 1st one is higher than 80%, only based on these profile will be built
	sequenceIdentity.push_back(true);		//to build the first profile based only on the first sequence
	std::string alNoLower;
	std::string alWithLower;
	for (int i = 1; i < seqNr; i++){
		alignPairwise(alNoLower,alWithLower,sequences.at(i).at(1),outputProfile,penalty, extensionPenalty,i,verbose);
		alignmentWithoutLowercase.push_back(alNoLower);
		alignmentWithLowercase.push_back(alWithLower);
		double identity = calcIdentity(alNoLower);
		if (identity > 0){
			sequenceIdentity.push_back(true);
		}
		else sequenceIdentity.push_back(false);	
		outputProfile.buildPseudoProfile(alignmentWithoutLowercase, sequenceIdentity);
	}
	return alignmentWithLowercase;		
}
//function performMSA for ENCODED SEQUENCES
std::vector<std::string> Sequences::performMSAencoded(Profile& outputProfile, FeaturesProfile& outputFeaturesProfile, double penalty, double extensionPenalty, std::string verbose, bool weightsModeOn, int dom, int phosph, int codon_length, std::vector<double>& identities){
	outputProfile = Profile(substitutionMatrix::convertToProfileFormat(sequencesEncoded.at(0).at(1))); 
	std::vector< std::vector<std::string> > alignmentWithoutLowercase;	//working alignment - without lowercase around cut out residues - would make latter aligning more complicated
	std::vector< std::vector<std::string> > alignmentWithLowercase;		//lowercase before and after cut out residues -- final result 
	alignmentWithoutLowercase.push_back(sequencesEncoded.at(0).at(1));
	alignmentWithLowercase.push_back(sequencesEncoded.at(0).at(1));
	std::vector<bool> sequenceIdentity; 					//'true' stored for every sequence which identity with the 1st one is higher than 80%, only based on these profile will be built
	identities.push_back(1); // identity of the 1st one to itself
	sequenceIdentity.push_back(true);					//to build the first profile based only on the first seqeunce
	outputFeaturesProfile.expandListOfFeatures(sequencesEncoded.at(0).at(1));
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase,sequenceIdentity,identities,weightsModeOn); 	//create features profile based on the 1st seq
	std::vector<std::string> alNoLower; //pairwise alignment without lowercase characters
	std::vector<std::string> alWithLower; //pairwise alignment with lowercase characters where chars were removed
	for (int i = 1; i < seqNr; i++){
		alignPairwise(alNoLower,alWithLower,sequencesEncoded.at(i).at(1),outputProfile,outputFeaturesProfile,penalty,extensionPenalty,i,verbose,codon_length);
		alignmentWithoutLowercase.push_back(alNoLower);
		alignmentWithLowercase.push_back(alWithLower);
		double identity = calcIdentity(alNoLower);
		identities.push_back(identity);
		if (identity > 0.95){
			sequenceIdentity.push_back(true);
		}
		else sequenceIdentity.push_back(false);	
		if (verbose!="0") outputProfile.printProfile();
	}
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase,sequenceIdentity,identities,weightsModeOn); 	//create features profile based on the 1st seq
	outputProfile.buildPseudoProfile(alignmentWithoutLowercase, sequenceIdentity,identities,weightsModeOn);
	return vecUtil::flatten(alignmentWithLowercase);		
}
//perform next round of MSA (good for all rounds except for the first one - you need a profile)
std::vector<std::string> Sequences::performMSAnextRound(Profile& outputProfile,FeaturesProfile& outputFeaturesProfile, double penalty, double extensionPenalty,std::string verbose, bool weightsModeOn, double identityCutoff,int dom, int phosph, int codon_length, std::vector<double> identities){
	std::vector< std::vector<std::string> > alignmentWithoutLowercase;	//working alignment - without lowercase around cut out residues - would make latter aligning more complicated
	std::vector< std::vector<std::string> > alignmentWithLowercase;		//lowercase before and after cut out residues -- final result 
	alignmentWithoutLowercase.push_back(sequencesEncoded.at(0).at(1));
	alignmentWithLowercase.push_back(sequencesEncoded.at(0).at(1));
	std::vector<std::string> alNoLower;
	std::vector<std::string> alWithLower;
	std::vector<bool> sequenceIdentity; 				//'true' stored for every sequence which identity with the 1st one is higher than identityCutoff, only based on these profile will be built
	sequenceIdentity.push_back(true);				//to build the first profile based only on the first seqeunce
	for (int i = 1; i < seqNr; i++){
		alignPairwise(alNoLower,alWithLower,sequencesEncoded.at(i).at(1),outputProfile,outputFeaturesProfile,penalty,extensionPenalty,i,verbose, codon_length);
		alignmentWithoutLowercase.push_back(alNoLower);
		alignmentWithLowercase.push_back(alWithLower);
		if (identities.at(i) > identityCutoff){
			sequenceIdentity.push_back(true);
		}
		else sequenceIdentity.push_back(false);	
	}
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase, sequenceIdentity,identities, weightsModeOn); 	//create features profile based on the 1st seq
	outputProfile.buildPseudoProfile(alignmentWithoutLowercase, sequenceIdentity, identities, weightsModeOn);
	if (verbose!="0") outputProfile.printProfile();
	return vecUtil::flatten(alignmentWithLowercase);		
}
//function calcIdentity
double Sequences::calcIdentity(const std::string& alignedSequence){
	double identicalResidues=0;
	for (int i = 0; i < alignedSequence.size(); i++){
		if (alignedSequence[i]==sequences.at(0)[1][i]){
			identicalResidues++;
		}
	}
	return identicalResidues/double(firstSequenceSize);
}
//function calcIdentity for ENCODED SEQUENCES
double Sequences::calcIdentity(const std::vector<std::string>& alignedSequence){
	double identicalResidues=0;
	for (int i = 0; i < alignedSequence.size(); i++){
		if (alignedSequence[i][0]==sequencesEncoded.at(0)[1][i][0]){
			identicalResidues++;
		}
	}
	return identicalResidues/double(firstSequenceSize);
}
//function getSequences
std::vector< std::vector<std::string> > Sequences::getSequences(){
	return sequences;
}
std::vector< std::vector<std::vector<std::string> > > Sequences::getEncodedSequences(){
	return sequencesEncoded;
}
//function removeGaps - takes pairwise alignment vector<string>, removes characters from the 2nd sequence that match gaps from 1st seq and returns vector<string> of 2 elements, where the 1st one is 2nd sequence with cut out chars and 2nd one is 2nd sequence with cut out chars and lowercase chars before and after that
void Sequences::removeGaps(std::string &alignmentWithLowercase, std::string &alignmentWithoutLowercase, std::vector<std::string>& alignment){
	std::vector<std::string> result;
	std::string s1 = alignment.at(0);
	std::string s2 = alignment.at(1);
	std::string newS2 = "";
	std::string newS2lower = "";
	char gap = '-';
	bool lowerFlag = false;
	for (int i = 0; i < alignment.at(0).size(); i++){
		char s1char = s1[i];
		if (s1char == gap){
			if (newS2lower.size() > 0){
				newS2lower.at(newS2lower.size()-1) = tolower(newS2lower.at(newS2lower.size()-1)); //change previous character to lowercase
			}
			lowerFlag = true; // flag to true so that the next character is also lowercase
		}
		else{
			if (lowerFlag){
				newS2lower += tolower(s2.at(i));
				newS2 += s2.at(i);
				lowerFlag = false;
			}
			else{
				newS2lower += s2.at(i);
				newS2 += s2.at(i);
			}
		}
	}
	alignmentWithLowercase = newS2lower;
	alignmentWithoutLowercase = newS2;
}
//function removeGaps - takes pairwise alignment vector<string>, removes characters from the 2nd sequence that match gaps from 1st seq and returns vector<string> of 2 elements, where the 1st one is 2nd sequence with cut out chars and 2nd one is 2nd sequence with cut out chars and lowercase chars before and after that ENCODED SEQUENCES
void Sequences::removeGaps(std::vector<std::string> &alignmentWithLowercase, std::vector<std::string> &alignmentWithoutLowercase, std::vector< std::vector<std::string> >& alignment){
	std::vector<std::vector<std::string> >result;
	std::vector<std::string> s1 = alignment.at(0);
	std::vector<std::string> s2 = alignment.at(1);
	std::vector<std::string> newS2;
	std::vector<std::string> newS2lower;
	char gap = '-';
	bool lowerFlag = false;
	for (int i = 0; i < alignment.at(0).size(); i++){
		char s1char = s1[i][0];
		if (s1char == gap){
			if (newS2lower.size() > 0){
				newS2lower.at(newS2lower.size()-1)[0] = tolower(newS2lower.at(newS2lower.size()-1)[0]); //change previous character to lowercase
			}
			lowerFlag = true; // flag to true so that the next character is also lowercase
		}
		else{
			if (lowerFlag){
				std::string newRes = s2.at(i);
				newRes[0] = tolower(newRes[0]);
				newS2lower.push_back(newRes);
				newS2.push_back(s2.at(i));
				lowerFlag = false;
			}
			else{
				newS2lower.push_back(s2.at(i));
				newS2.push_back(s2.at(i));
			}
		}
	}
	alignmentWithLowercase = newS2lower;
	alignmentWithoutLowercase = newS2;
}
//function alignPairwise
void Sequences::alignPairwise(std::string &alNoLower,std::string &alWithLower,std::string seq2, Profile& prf, double penalty, double extensionPenalty, int deb, std::string verbose){
	int profileLength = prf.getMatrix().at(0).size();
	std::vector<std::string> alignment;
	ScoringMatrix scores(profileLength,seq2.size(),penalty, extensionPenalty);
	scores.calculateScores(seq2, prf, deb);
	scores.nwAlignment(&alignment,seq2,prf,verbose);
	removeGaps(alWithLower,alNoLower,alignment); 
}
//function alignPairwise for ENCODED SEQUENCES
void Sequences::alignPairwise(std::vector<std::string> &alNoLower,std::vector<std::string> &alWithLower,std::vector<std::string> seq2, Profile& prf, FeaturesProfile& featPrf, double penalty, double extensionPenalty,int deb, std::string verbose, int codon_length){
	int profileLength = prf.getMatrix().at(0).size();
	std::vector< std::vector<std::string> > alignment;
	ScoringMatrix scores(profileLength,seq2.size(),penalty,extensionPenalty);
	scores.calculateScores(seq2, prf, featPrf,deb);
	scores.nwAlignment(&alignment,seq2,prf,featPrf,verbose, codon_length);
	removeGaps(alWithLower,alNoLower,alignment); 
}
