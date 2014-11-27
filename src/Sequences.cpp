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

Sequences::Sequences(std::vector<std::vector<std::vector<std::string> > >& s){
	std::vector<std::string> additional_features;
	for (unsigned int i = 0; i < s.size(); i++){
		//sequence_names.push_back(s[i][0]);
		m_sequence_names.push_back(s[i][0][0]);
		sequence new_seq;
		for (unsigned int j = 0; j < s[i][1].size(); j++){
			Residue newRes(s[i][1][j], additional_features);
			new_seq.push_back(newRes);
		}
		m_sequences_aa.push_back(new_seq);
	}
	m_seqNr = s.size(); 
	m_firstSequenceSize = s[0][1].size();
}


Sequences::Sequences(){
}


//function performMSAfirstround - performs the first round of alignments, 
//all vs query seq (first calculates profile based only on the query seq, then 
//aligns all sequences and calculates identity of each sequence to the query seq.)
std::vector<std::string> Sequences::performMSAfirstround(Profile& outputProfile, 
                                                         FeaturesProfile& outputFeaturesProfile, 
                                                         double penalty, 
                                                         double endPenalty, 
                                                         double extensionPenalty, 
                                                         bool weightsModeOn, 
                                                         int codon_length, 
                                                         std::vector<double>& identities){
	outputProfile = Profile(substitutionMatrix::convertToProfileFormat(m_sequences_aa[0])); 
  //working alignment - without lowercase around cut out residues
	sequenceList alignmentWithoutLowercase;	
  //lowercase before and after cut out residues -- final result 
	sequenceList alignmentWithLowercase;		
	alignmentWithoutLowercase.push_back(m_sequences_aa[0]);
	alignmentWithLowercase.push_back(m_sequences_aa[0]);
  //'true' stored for every sequence which identity with the 1st one is higher 
  //than 80%, only based on these profile will be built
	std::vector<bool> sequenceIdentity; 					
  // identity of the 1st one to itself
	identities.push_back(1); 
  //to build the first profile based only on the first seqeunce
	sequenceIdentity.push_back(true);					
	outputFeaturesProfile.expandListOfFeatures(m_sequences_aa);
  //create features profile based on the 1st seq
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase,
                                      identities, weightsModeOn, codon_length); 	
  add_feature_indexes(outputFeaturesProfile);
  //pairwise alignment without lowercase characters
	sequence alNoLower; 
  //pairwise alignment with lowercase characters where chars were removed
	sequence alWithLower; 
  for (auto &seqI: m_sequences_aa){
		alignPairwise(alNoLower, alWithLower, seqI, outputProfile, 
                  outputFeaturesProfile, penalty, endPenalty, extensionPenalty, 
                  0, codon_length);
		double identity = calcIdentity(alNoLower);
		identities.push_back(identity);
		if (identity > 0.9){
			alignmentWithoutLowercase.push_back(alNoLower);
			alignmentWithLowercase.push_back(alWithLower);
		}
		else sequenceIdentity.push_back(false);	
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
                                     bool weightsModeOn, 
                                     double identityCutoff,
                                     int codon_length, 
                                     std::vector<double>& identities, 
                                     int& prev_alignments){
	int next_alignments = countAlignments(identityCutoff, identities);
	if (next_alignments > prev_alignments){
    //working alignment - without lowercase around cut out residues
    //would make latter aligning more complicated
		sequenceList alignmentWithoutLowercase;	
    //lowercase before and after cut out residues -- final result 
		sequenceList alignmentWithLowercase;		
		alignmentWithoutLowercase.push_back(m_sequences_aa[0]);
		alignmentWithLowercase.push_back(m_sequences_aa[0]);
    // tmp pairwise alignment (and so is alWithLower)
		sequence alNoLower; 
		sequence alWithLower;
		for (int i = 1; i < m_seqNr; i++){
			if (identities[i] > identityCutoff){
        // NW alignment of the ith seq against the profile
				alignPairwise(alNoLower, alWithLower, m_sequences_aa[i], outputProfile,
                      outputFeaturesProfile,penalty,endPenalty,extensionPenalty,
                      0, codon_length); 
				alignmentWithoutLowercase.push_back(alNoLower);
				alignmentWithLowercase.push_back(alWithLower);
			}
		}
    //create features profile based on the 1st seq
		outputFeaturesProfile.createProfile(alignmentWithoutLowercase,identities, 
                                        weightsModeOn, codon_length);
		outputProfile.buildPseudoProfile(alignmentWithoutLowercase,identities,
                                     weightsModeOn);
		*prevAlignment = vecUtil::flatten(alignmentWithLowercase);
		prev_alignments = next_alignments; //number of performed alignments
	}
}


//function calcIdentity, calculates identity with the query sequence (takes 
//aligned sequence with the gaps cut out)
double Sequences::calcIdentity(const sequence& alignedSequence){
	double identicalResidues=0;
	for (unsigned int i = 0; i < alignedSequence.size(); i++){
		if (alignedSequence[i].getAA() == m_sequences_aa[0][i].getAA()){
			identicalResidues++;
		}
	}
	return identicalResidues/double(m_firstSequenceSize);
}


//function removeGaps - takes pairwise alignment vector<string>, removes 
//characters from the 2nd sequence that match gaps from 1st seq and returns 
//vector<string> of 2 elements, where the 1st one is 2nd sequence with cut out 
//chars and 2nd one is 2nd sequence with cut out chars and lowercase chars 
//before and after that
void Sequences::removeGaps(sequence& alignmentWithLowercase, 
                           sequence& alignmentWithoutLowercase, 
                           sequenceList& alignment){
	sequence s1 = alignment[0];
	sequence s2 = alignment[1];
	sequence newS2;
	sequence newS2lower;
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
void Sequences::alignPairwise(sequence& alNoLower, 
                              sequence& alWithLower, 
                              sequence& seq2, 
                              Profile& prf, 
                              FeaturesProfile& featPrf, 
                              double penalty, double endPenalty, 
                              double extensionPenalty, int deb, 
                              int codon_length){
	int profileLength = prf.getMatrix()[0].size();
	sequenceList alignment;
	ScoringMatrix scores(profileLength, seq2.size(), penalty, 
                       endPenalty, extensionPenalty);
	scores.calculateScores(seq2, prf, featPrf, deb, 
                         codon_length);
	scores.nwAlignment(&alignment, seq2, prf, featPrf, 
                     codon_length);

	removeGaps(alWithLower,alNoLower,alignment); 
}


//count alignments that will be performed in this round
int Sequences::countAlignments(double identity_cutoff, 
                               std::vector<double>& identities){
	int count = 0;
	//for (unsigned int i = 0; i < identities.size(); i++){
  for (auto &item: identities){
		if (item > identity_cutoff){
			count++;
		}
	}
	return count;
}


//adds features from the tuple 'feature_rules'(usr defined) to relevant 
//residues (also specified in 'feature_rules')
void Sequences::add_usr_features(std::vector<std::tuple<std::string,std::string, 
                                 int, int, int, double, double, double, double, 
                                 std::string, std::string> >& feature_rules){
	//for (unsigned int i = 0; i < feature_rules.size(); i++){
	for (auto &rule: feature_rules){
		std::string feat_name = std::string("USR_")
                            + std::get<0>(rule)
                            + std::string("_")
                            + std::get<1>(rule);
		int sequence_no = std::get<2>(rule);
		int start = std::get<3>(rule);
		int end = std::get<4>(rule)+1;
    signed int seq_length = m_sequences_aa[sequence_no].size();
    if (sequence_no < (signed)m_sequences_aa.size()){
		  for (int j = start; j < end && j < seq_length; j++){
		  	m_sequences_aa[sequence_no][j].add_feature(feat_name);
		  }
    }
	}
}


std::vector< std::string> Sequences::get_names(){
  return m_sequence_names;
}

void Sequences::add_feature_indexes(FeaturesProfile& fprf){
  std::string nothing = "AA";
  for (auto &seq: m_sequences_aa){
    for (auto &res: seq){
        std::vector<std::string> features = res.getFeatures();
        std::vector<int> indexes;
        for (auto &feat: features){
          if (feat != nothing){
            indexes.push_back(fprf.findFeaturesIndex(feat));
          }
        }
        res.setFeatIndexes(indexes);
    }
  }
}
