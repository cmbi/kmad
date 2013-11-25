//Sequences class implementation
#include <iostream>
#include <string>
#include <vector>
#include "Sequences.h"
#include "Profile.h"
#include "ScoringMatrix.h"
#include "UsefulStuff.h"
#include "SubstitutionMatrix.h"
using namespace std;
extern UsefulStuff util;
// 1st constructor - assigns s to vector< vector<string> > sequences - [[seqName,sequence],...]
Sequences::Sequences(vector< vector<string> > s){
	sequences = s;
	seqNr = s.size(); 
	firstSequenceSize = s.at(0).at(1).size();
}
Sequences::Sequences(vector< vector< vector<string> > > s){
	sequencesEncoded = s;
	seqNr = s.size(); 
	firstSequenceSize = s.at(0).at(1).size();
}
Sequences::Sequences(){
}
//function performMSA - returns ready alignment, with gaps cut out from the first sequence and so on...
vector<string> Sequences::performMSA(Profile *outputProfile,int penalty, bool verbose){
	SubstitutionMatrix blosum = SubstitutionMatrix();
	Profile pseudoProfile(blosum.convertToProfileFormat(sequences.at(0).at(1))); //IMPLEMENT convertToProfileFormat // psuedoProfile - it's the one that will later be combined profile and sbst matrix
	vector<string> alignmentWithoutLowercase;	//working alignment - without lowercase around cut out residues - would make latter aligning more complicated
	vector<string> alignmentWithLowercase;		//lowercase before and after cut out residues -- final result 
	alignmentWithoutLowercase.push_back(sequences.at(0).at(1));
	alignmentWithLowercase.push_back(sequences.at(0).at(1));
	vector<string> tmpAlignment;
	vector<bool> sequenceIdentity; 			//'true' stored for every sequence which identity with the 1st one is higher than 80%, only based on these profile will be built
	sequenceIdentity.push_back(true);		//we always want to build profile based on the first seqeunce
	string alNoLower;
	string alWithLower;
	bool flag = false;
	for (int i = 1; i < seqNr; i++){
		alignPairwise(&alNoLower,&alWithLower,sequences.at(i).at(1),pseudoProfile,penalty,i,verbose);
		alignmentWithoutLowercase.push_back(alNoLower);
		alignmentWithLowercase.push_back(alWithLower);
		double identity = calcIdentity(alNoLower);
		if (identity > 0.5){
			sequenceIdentity.push_back(true);
		}
		else sequenceIdentity.push_back(false);	
		pseudoProfile.buildPseudoProfile(alignmentWithoutLowercase, sequenceIdentity, blosum);
		if (verbose) pseudoProfile.printProfile();
	}
	*outputProfile = pseudoProfile;
	return alignmentWithLowercase;		
}
//function performMSAencoded
vector<string> Sequences::performMSAencoded(Profile *outputProfile,Profile *outputFeaturesProfile, int penalty, bool verbose){
	SubstitutionMatrix blosum = SubstitutionMatrix();
	Profile pseudoProfile(blosum.convertToProfileFormat(sequences.at(0).at(1))); //IMPLEMENT convertToProfileFormat // psuedoProfile - it's the one that will later be combined profile and sbst matrix
	vector<string> alignmentWithoutLowercase;	//working alignment - without lowercase around cut out residues - would make latter aligning more complicated
	vector<string> alignmentWithLowercase;		//lowercase before and after cut out residues -- final result 
	alignmentWithoutLowercase.push_back(sequences.at(0).at(1));
	alignmentWithLowercase.push_back(sequences.at(0).at(1));
	vector<string> tmpAlignment;
	vector<bool> sequenceIdentity; 			//'true' stored for every sequence which identity with the 1st one is higher than 80%, only based on these profile will be built
	sequenceIdentity.push_back(true);		//we always want to build profile based on the first seqeunce
	string alNoLower;
	string alWithLower;
	bool flag = false;
	for (int i = 1; i < seqNr; i++){
		alignPairwise(&alNoLower,&alWithLower,sequences.at(i).at(1),pseudoProfile,penalty,i,verbose);
		alignmentWithoutLowercase.push_back(alNoLower);
		alignmentWithLowercase.push_back(alWithLower);
		double identity = calcIdentity(alNoLower);
		if (identity > 0.5){
			sequenceIdentity.push_back(true);
		}
		else sequenceIdentity.push_back(false);	
		pseudoProfile.buildPseudoProfile(alignmentWithoutLowercase, sequenceIdentity, blosum);
		if (verbose) pseudoProfile.printProfile();
	}
	*outputProfile = pseudoProfile;
	return alignmentWithLowercase;		
}
//function calcIdentity
double Sequences::calcIdentity(const string& alignedSequence){
	double identicalResidues=0;
	for (int i = 0; i < alignedSequence.size(); i++){
		if (alignedSequence[i]==sequences.at(0)[1][i]){
			identicalResidues++;
		}
	}
	return identicalResidues/double(firstSequenceSize);
}
//function getSequences
vector< vector<string> > Sequences::getSequences(){
	return sequences;
}
//function removeGaps - takes pairwise alignment vector<string>, removes characters from the 2nd sequence that match gaps from 1st seq and returns vector<string> of 2 elements, where the 1st one is 2nd sequence with cut out chars and 2nd one is 2nd sequence with cut out chars and lowercase chars before and after that
void Sequences::removeGaps(string *alignmentWithLowercase, string *alignmentWithoutLowercase, vector<string>& alignment){
	vector<string> result;
	string s1 = alignment.at(0);
	string s2 = alignment.at(1);
	string newS2 = "";
	string newS2lower = "";
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
	*alignmentWithLowercase = newS2lower;
	*alignmentWithoutLowercase = newS2;
}
//function alignPairwise
void Sequences::alignPairwise(string *alNoLower,string *alWithLower,string seq2, Profile& prf, int penalty, int deb, bool verbose){
	int profileLength = prf.getMatrix().at(0).size();
	string tmpResultWithLowercase;
	string tmpResultWithoutLowercase;
	vector<string> alignment;
	ScoringMatrix scores(profileLength,seq2.size(),penalty);
	scores.calculateScores(seq2, prf, deb);
	scores.nwAlignment(&alignment,seq2,prf,verbose);
	removeGaps(&tmpResultWithLowercase,&tmpResultWithoutLowercase,alignment); 
	*alNoLower = tmpResultWithoutLowercase;
	*alWithLower = tmpResultWithLowercase;
}
void Sequences::printEncodedSequence(int seqNo){
	cout << sequencesEncoded.at(seqNo).at(0).at(0)<<"\n";
	vector<string> tmpSequence = sequencesEncoded.at(seqNo).at(1);
	for (int i = 0; i < tmpSequence.size();i++){
		cout << tmpSequence.at(i);
	}
	cout << "\n";
}
