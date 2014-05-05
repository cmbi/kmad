#include "Sequences.h"
#include "FeaturesProfile.h"
#include "Profile.h"
#include "ScoringMatrix.h"
#include "substitutionMatrix.h"
#include "vecUtil.h"
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
//constructor for ENCODED SEQUENCES
Sequences::Sequences(std::vector< std::vector< std::vector<std::string> > > s){
	sequencesEncoded = s;
	seqNr = s.size(); 
	firstSequenceSize = s.at(0).at(1).size();
}
//function performMSA for ENCODED SEQUENCES
std::vector<std::string> Sequences::performMSAencoded(Profile& outputProfile, FeaturesProfile& outputFeaturesProfile, double penalty, double extensionPenalty, std::string verbose, bool weightsModeOn, int codon_length, std::vector<double>& identities){
	outputProfile = Profile(substitutionMatrix::convertToProfileFormat(sequencesEncoded.at(0).at(1))); 
	std::vector< std::vector<std::string> > alignmentWithoutLowercase;	//working alignment - without lowercase around cut out residues - would make latter aligning more complicated
	std::vector< std::vector<std::string> > alignmentWithLowercase;		//lowercase before and after cut out residues -- final result 
	alignmentWithoutLowercase.push_back(sequencesEncoded.at(0).at(1));

	alignmentWithLowercase.push_back(sequencesEncoded.at(0).at(1));
	std::vector<bool> sequenceIdentity; 					//'true' stored for every sequence which identity with the 1st one is higher than 80%, only based on these profile will be built
	identities.push_back(1); // identity of the 1st one to itself
	sequenceIdentity.push_back(true);					//to build the first profile based only on the first seqeunce
	outputFeaturesProfile.expandListOfFeatures(sequencesEncoded.at(0).at(1), codon_length);
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase,identities,weightsModeOn,codon_length); 	//create features profile based on the 1st seq

	std::vector<std::string> alNoLower; //pairwise alignment without lowercase characters
	std::vector<std::string> alWithLower; //pairwise alignment with lowercase characters where chars were removed
	for (int i = 1; i < seqNr; i++){
		time_t start = clock(); 
		alignPairwise(alNoLower,alWithLower,sequencesEncoded.at(i).at(1),outputProfile,outputFeaturesProfile,penalty,extensionPenalty,i,verbose,codon_length);
		time_t end = clock(); 
		//std::cout << "Time clock() = " << (end - start)/(double)CLOCKS_PER_SEC << std::endl;
		double identity = calcIdentity(alNoLower);
		identities.push_back(identity);
		if (identity > 0.9){
			alignmentWithoutLowercase.push_back(alNoLower);
			alignmentWithLowercase.push_back(alWithLower);
		}
		else sequenceIdentity.push_back(false);	
		if (verbose!="0") outputProfile.printProfile();
	}
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase,identities,weightsModeOn,codon_length); 	//create features profile based on the 1st seq
	outputProfile.buildPseudoProfile(alignmentWithoutLowercase,identities,weightsModeOn);
	return vecUtil::flatten(alignmentWithLowercase);		
}
//perform next round of MSA (good for all rounds except for the first one - you need a profile)
std::vector<std::string> Sequences::performMSAnextRound(Profile& outputProfile,FeaturesProfile& outputFeaturesProfile, double penalty, double extensionPenalty,std::string verbose, bool weightsModeOn, double identityCutoff,int codon_length, std::vector<double> identities){
	std::vector< std::vector<std::string> > alignmentWithoutLowercase;	//working alignment - without lowercase around cut out residues - would make latter aligning more complicated
	std::vector< std::vector<std::string> > alignmentWithLowercase;		//lowercase before and after cut out residues -- final result 
	alignmentWithoutLowercase.push_back(sequencesEncoded.at(0).at(1));
	alignmentWithLowercase.push_back(sequencesEncoded.at(0).at(1));
	std::vector<std::string> alNoLower;
	std::vector<std::string> alWithLower;
	std::vector<bool> sequenceIdentity; 				//'true' stored for every sequence which identity with the 1st one is higher than identityCutoff, only based on these profile will be built
	sequenceIdentity.push_back(true);				//to build the first profile based only on the first seqeunce
	int iter=0;
	for (int i = 1; i < seqNr; i++){
		if (identities.at(i) > identityCutoff){
			alignPairwise(alNoLower,alWithLower,sequencesEncoded.at(i).at(1),outputProfile,outputFeaturesProfile,penalty,extensionPenalty,i,verbose, codon_length);
			alignmentWithoutLowercase.push_back(alNoLower);
			alignmentWithLowercase.push_back(alWithLower);
		}
	}
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase,identities, weightsModeOn, codon_length);//create features profile based on the 1st seq
	outputProfile.buildPseudoProfile(alignmentWithoutLowercase,identities, weightsModeOn);
	if (verbose!="0") outputProfile.printProfile();
	return vecUtil::flatten(alignmentWithLowercase);		
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
//function alignPairwise for ENCODED SEQUENCES
void Sequences::alignPairwise(std::vector<std::string> &alNoLower,std::vector<std::string> &alWithLower,std::vector<std::string> seq2, Profile& prf, FeaturesProfile& featPrf, double penalty, double extensionPenalty,int deb, std::string verbose, int codon_length){
	int profileLength = prf.getMatrix().at(0).size();
	std::vector< std::vector<std::string> > alignment;
	time_t p1 = clock();
	ScoringMatrix scores(profileLength,seq2.size(),penalty,extensionPenalty);
	time_t p2 = clock();
	scores.calculateScores(seq2, prf, featPrf,deb, codon_length);
	time_t p3 = clock();
	scores.nwAlignment(&alignment,seq2,prf,featPrf,verbose, codon_length);
	time_t p4 = clock();
	removeGaps(alWithLower,alNoLower,alignment); 
	time_t p5 = clock();
	/*
	std::cout << "constructor = " << (p2 - p1)/(double)CLOCKS_PER_SEC << std::endl;
	std::cout << "calcScores = " << (p3 - p2)/(double)CLOCKS_PER_SEC << std::endl;
	std::cout << "nwALignment = " << (p4 - p3)/(double)CLOCKS_PER_SEC << std::endl;
	std::cout << "removeGaps = " << (p5 - p4)/(double)CLOCKS_PER_SEC << std::endl;
	*/
}
