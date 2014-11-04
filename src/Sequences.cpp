#include "Sequences.h"
#include "Residue.h"
#include "FeaturesProfile.h"
#include "Profile.h"
#include "ScoringMatrix.h"
#include "substitutionMatrix.h"
#include "vecUtil.h"
#include <iostream>
#include <tuple>
#include <string>
#include <vector>
#include <ctime>
//constructor for ENCODED SEQUENCES
Sequences::Sequences(std::vector< std::vector< std::vector<std::string> > > s){
	sequencesEncoded = s;
	std::vector<std::string> additional_features;
	for (unsigned int i = 0; i < sequencesEncoded.size(); i++){
		sequence_names.push_back(sequencesEncoded[i][0]);
		std::vector<Residue> new_seq;
		for (unsigned int j = 0; j < sequencesEncoded[i][1].size(); j++){
			Residue newRes(sequencesEncoded[i][1][j], additional_features);
			new_seq.push_back(newRes);
		}
		sequences_aa.push_back(new_seq);
	}
	seqNr = s.size(); 
	firstSequenceSize = s[0][1].size();
}
//function performMSAfirstround - performs the first round of alignments, 
//all vs query seq (first calculates profile based only on the query seq, then 
//aligns all sequences and calculates identity of each sequence to the query seq.)
std::vector<std::string> Sequences::performMSAfirstround(Profile& outputProfile, 
                                                         FeaturesProfile& outputFeaturesProfile, 
                                                         double penalty, 
                                                         double endPenalty, 
                                                         double extensionPenalty, 
                                                         std::string verbose, 
                                                         bool weightsModeOn, 
                                                         int codon_length, 
                                                         std::vector<double>& identities){
	outputProfile = Profile(substitutionMatrix::convertToProfileFormat(sequences_aa[0])); 
  //working alignment - without lowercase around cut out residues
	std::vector< std::vector<Residue> > alignmentWithoutLowercase;	
  //lowercase before and after cut out residues -- final result 
	std::vector< std::vector<Residue> > alignmentWithLowercase;		
	alignmentWithoutLowercase.push_back(sequences_aa[0]);
	alignmentWithLowercase.push_back(sequences_aa[0]);
  //'true' stored for every sequence which identity with the 1st one is higher 
  //than 80%, only based on these profile will be built
	std::vector<bool> sequenceIdentity; 					
  // identity of the 1st one to itself
	identities.push_back(1); 
  //to build the first profile based only on the first seqeunce
	sequenceIdentity.push_back(true);					
	outputFeaturesProfile.expandListOfFeatures(sequences_aa);
  //create features profile based on the 1st seq
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase,
                                      identities, weightsModeOn, codon_length); 	
  //pairwise alignment without lowercase characters
	std::vector<Residue> alNoLower; 
  //pairwise alignment with lowercase characters where chars were removed
	std::vector<Residue> alWithLower; 
	for (int i = 1; i < seqNr; i++){
		alignPairwise(alNoLower, alWithLower, sequences_aa[i], outputProfile, 
                  outputFeaturesProfile, penalty, endPenalty, extensionPenalty, 
                  i, verbose, codon_length, i);
		double identity = calcIdentity(alNoLower);
		identities.push_back(identity);
		if (identity > 0.9){
			alignmentWithoutLowercase.push_back(alNoLower);
			alignmentWithLowercase.push_back(alWithLower);
		}
		else sequenceIdentity.push_back(false);	
		if (verbose!="0") outputProfile.printProfile();
	}
  //create features profile based on the 1st seq
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase, 
                                      identities, 
                                      weightsModeOn, 
                                      codon_length);
	outputProfile.buildPseudoProfile(alignmentWithoutLowercase,
                                   identities,
                                   weightsModeOn);
	return vecUtil::flatten(alignmentWithLowercase);		
}
//perform next round of MSA (good for all rounds except for the first one - 
//you need a profile)
void Sequences::performMSAnextRounds(std::vector<std::string>* prevAlignment, 
                                     Profile& outputProfile,
                                     FeaturesProfile& outputFeaturesProfile, 
                                     double penalty, 
                                     double endPenalty, 
                                     double extensionPenalty,
                                     std::string verbose, 
                                     bool weightsModeOn, 
                                     double identityCutoff,
                                     int codon_length, 
                                     std::vector<double>& identities, 
                                     int& prev_alignments){
	int next_alignments = countAlignments(identityCutoff, identities);

	if (next_alignments > prev_alignments){
    //working alignment - without lowercase around cut out residues
    //would make latter aligning more complicated
		std::vector< std::vector<Residue> > alignmentWithoutLowercase;	
    //lowercase before and after cut out residues -- final result 
		std::vector< std::vector<Residue> > alignmentWithLowercase;		
		alignmentWithoutLowercase.push_back(sequences_aa[0]);
		alignmentWithLowercase.push_back(sequences_aa[0]);
    // tmp pairwise alignment (and so is alWithLower)
		std::vector<Residue> alNoLower; 
		std::vector<Residue> alWithLower;
		for (int i = 1; i < seqNr; i++){
			if (identities[i] > identityCutoff){
        // NW alignment of the ith seq against the profile
				alignPairwise(alNoLower, alWithLower, sequences_aa[i], outputProfile,
                      outputFeaturesProfile,penalty,endPenalty,extensionPenalty,
                      i, verbose, codon_length, i); 
				alignmentWithoutLowercase.push_back(alNoLower);
				alignmentWithLowercase.push_back(alWithLower);
			}
		}
    //create features profile based on the 1st seq
		outputFeaturesProfile.createProfile(alignmentWithoutLowercase,identities, 
                                        weightsModeOn, codon_length);
		outputProfile.buildPseudoProfile(alignmentWithoutLowercase,identities,
                                     weightsModeOn);
		if (verbose!="0") outputProfile.printProfile();
		*prevAlignment = vecUtil::flatten(alignmentWithLowercase);
		prev_alignments = next_alignments; //number of performed alignments
	}
}
//function calcIdentity, calculates identity with the query sequence (takes 
//aligned sequence with the gaps cut out)
double Sequences::calcIdentity(const std::vector<Residue>& alignedSequence){
	double identicalResidues=0;
	for (unsigned int i = 0; i < alignedSequence.size(); i++){
		if (alignedSequence[i].getAA() == sequences_aa[0][i].getAA()){
			identicalResidues++;
		}
	}
	return identicalResidues/double(firstSequenceSize);
}
std::vector< std::vector<std::vector<std::string> > > Sequences::getEncodedSequences(){
	return sequencesEncoded;
}
//function removeGaps - takes pairwise alignment vector<string>, removes 
//characters from the 2nd sequence that match gaps from 1st seq and returns 
//vector<string> of 2 elements, where the 1st one is 2nd sequence with cut out 
//chars and 2nd one is 2nd sequence with cut out chars and lowercase chars 
//before and after that
void Sequences::removeGaps(std::vector<Residue> &alignmentWithLowercase, 
                           std::vector<Residue> &alignmentWithoutLowercase, 
                           std::vector< std::vector<Residue> >& alignment){
	std::vector<std::vector<Residue> >result;
	std::vector<Residue> s1 = alignment[0];
	std::vector<Residue> s2 = alignment[1];
	std::vector<Residue> newS2;
	std::vector<Residue> newS2lower;
	char gap = '-';
	bool lowerFlag = false;
	for (unsigned int i = 0; i < alignment[0].size(); i++){
		char s1char = s1[i].getAA();
		if (s1char == gap){
			if (newS2lower.size() > 0){
        //change previous character to lowercase
				newS2lower[newS2lower.size()-1].lowercase(); 
			}
      // flag to true so that the next character is also lowercase
			lowerFlag = true; 
		}
		else{
			if (lowerFlag){   //lowercase char
				Residue newRes = s2[i];
				newRes.lowercase();
        //add lowercase char to the alignment with lowercases
				newS2lower.push_back(newRes); 
        //add uppercase alignment to the alignment without lowercases
				newS2.push_back(s2[i]);      	
				lowerFlag = false;
			}
			else{ 		
        //uppercase char
        // adds the same uppercase char to both alignments (with lowercases and 
        // without lowercases)
				newS2lower.push_back(s2[i]); 
				newS2.push_back(s2[i]);
			}
		}
	}
	alignmentWithLowercase = newS2lower;
	alignmentWithoutLowercase = newS2;
}
//function alignPairwise -> takes a sequence and profiles, returns a ready 
//alignment of the two, with gaps cut out
void Sequences::alignPairwise(std::vector<Residue> &alNoLower, 
                              std::vector<Residue> &alWithLower, 
                              std::vector<Residue>& seq2, 
                              Profile& prf, 
                              FeaturesProfile& featPrf, 
                              double penalty, double endPenalty, 
                              double extensionPenalty, int deb, 
                              std::string verbose, int codon_length, 
                              int sequence_no){
	int profileLength = prf.getMatrix()[0].size();
	std::vector< std::vector<Residue> > alignment;
	ScoringMatrix scores(profileLength, seq2.size(), penalty, 
                       endPenalty, extensionPenalty);
	scores.calculateScores(seq2, prf, featPrf, deb, 
                         codon_length, sequence_no);
	scores.nwAlignment(&alignment, seq2, prf, featPrf, 
                     verbose, codon_length, sequence_no);

	removeGaps(alWithLower,alNoLower,alignment); 
}
//count alignments that will be performed in this round
int Sequences::countAlignments(double identity_cutoff, 
                               std::vector<double>& identities){
	int count = 0;
	for (unsigned int i = 0; i < identities.size(); i++){
		if (identities[i] > identity_cutoff){
			count++;
		}
	}
	return count;
}
void Sequences::printSequence(int seq_index) const{
	for (unsigned int i = 0; i < sequences_aa[seq_index].size(); i++){
		std::cout << sequences_aa[seq_index][i].getAA();
	}
	std::cout << std::endl;
}
//adds features from the tuple 'feature_rules'(usr defined) to relevant 
//residues (also specified in 'feature_rules')
void Sequences::add_usr_features(std::vector<std::tuple<std::string,std::string, 
                                 int, int, int, double, double, double, double, 
                                 std::string, std::string> >& feature_rules){
	for (unsigned int i = 0; i < feature_rules.size(); i++){
		std::string feat_name = std::string("USR_")
                            + std::get<0>(feature_rules[i])
                            + std::string("_")
                            + std::get<1>(feature_rules[i]);
		int sequence_no = std::get<2>(feature_rules[i]);
		int start = std::get<3>(feature_rules[i]);
		int end = std::get<4>(feature_rules[i])+1;
		for (int j = start; j < end; j++){
			sequences_aa[sequence_no][j].add_feature(feat_name);
		}
	}
}
