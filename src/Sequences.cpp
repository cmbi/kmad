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


Sequences::Sequences(codonSeqWithNamesList &s){
	for (unsigned int i = 0; i < s.size(); i++){
		m_sequence_names.push_back(s[i][0][0]);
		sequence new_seq;
		for (unsigned int j = 0; j < s[i][1].size(); j++){
			Residue newRes(s[i][1][j]);
			new_seq.push_back(newRes);
		}
		m_sequences_aa.push_back(new_seq);
	}
	m_seqNr = s.size(); 
	m_firstSequenceSize = s[0][1].size();
}


Sequences::Sequences(){
}


string_sequences Sequences::performMSAfirstround(Profile& outputProfile, 
                                                 FeaturesProfile& outputFeaturesProfile, 
                                                 double penalty, 
                                                 double endPenalty, 
                                                 double extensionPenalty, 
                                                 int codon_length, 
                                                 identitiesList& identities){
	outputProfile = Profile(substitutionMatrix::convertToProfileFormat(m_sequences_aa[0])); 
  //working alignment - without lowercase around cut out residues
	sequenceList alignmentWithoutLowercase;	
  //lowercase before and after cut out residues -- final result 
	sequenceList alignmentWithLowercase;		
	alignmentWithoutLowercase.push_back(m_sequences_aa[0]);
	alignmentWithLowercase.push_back(m_sequences_aa[0]);
  // identity of the 1st one to itself
  //to build the first profile based only on the first seqeunce
	identities.push_back(1); 
	outputFeaturesProfile.expandListOfFeatures(m_sequences_aa);
  //create features profile based on the 1st seq
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase,
                                      identities, codon_length); 	
  add_feature_indexes(outputFeaturesProfile);
  //pairwise alignment without lowercase characters
	sequence alNoLower; 
  //pairwise alignment with lowercase characters where chars were removed
	sequence alWithLower; 
  for (auto &seqI: m_sequences_aa){
		alignPairwise(alNoLower, alWithLower, seqI, outputProfile, 
                  outputFeaturesProfile, penalty, endPenalty, extensionPenalty, 
                  codon_length);
		double identity = calcIdentity(alNoLower);
		identities.push_back(identity);
		if (identity > 0.9){
			alignmentWithoutLowercase.push_back(alNoLower);
			alignmentWithLowercase.push_back(alWithLower);
		}
	}
  //create features profile based on the 1st seq
	outputFeaturesProfile.createProfile(alignmentWithoutLowercase, 
                                      identities, 
                                      codon_length);
	outputProfile.processProfile(alignmentWithoutLowercase,
                               identities);
	return vecUtil::flatten(alignmentWithLowercase);		
}


void Sequences::performMSAnextRounds(string_sequences& prevAlignment, 
                                     Profile& outputProfile,
                                     FeaturesProfile& outputFeaturesProfile, 
                                     double penalty, 
                                     double endPenalty, 
                                     double extensionPenalty,
                                     double identityCutoff,
                                     int codon_length, 
                                     identitiesList& identities, 
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
                      codon_length); 
				alignmentWithoutLowercase.push_back(alNoLower);
				alignmentWithLowercase.push_back(alWithLower);
			}
		}
    //create features profile based on the 1st seq
		outputFeaturesProfile.createProfile(alignmentWithoutLowercase,identities, 
                                        codon_length);
		outputProfile.processProfile(alignmentWithoutLowercase,identities);
		prevAlignment = vecUtil::flatten(alignmentWithLowercase);
    //update number of performed alignments
		prev_alignments = next_alignments; 
	}
}


double Sequences::calcIdentity(const sequence& alignedSequence){
	double identicalResidues=0;
	for (unsigned int i = 0; i < alignedSequence.size(); i++){
		if (alignedSequence[i].getAA() == m_sequences_aa[0][i].getAA()){
			identicalResidues++;
		}
	}
	return identicalResidues/double(m_firstSequenceSize);
}


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


void Sequences::alignPairwise(sequence& alNoLower, 
                              sequence& alWithLower, 
                              sequence& seq2, 
                              Profile& prf, 
                              FeaturesProfile& featPrf, 
                              double penalty, double endPenalty, 
                              double extensionPenalty, 
                              int codon_length){
	int profileLength = prf.getMatrix()[0].size();
	sequenceList alignment;
	ScoringMatrix scores(profileLength, seq2.size(), penalty, 
                       endPenalty, extensionPenalty);
	scores.calculateScores(seq2, prf, featPrf, 
                         codon_length);
	scores.nwAlignment(&alignment, seq2, prf, featPrf, 
                     codon_length);

	removeGaps(alWithLower,alNoLower,alignment); 
}


int Sequences::countAlignments(double identity_cutoff, 
                               identitiesList& identities){
	int count = 0;
  for (auto &item: identities){
		if (item > identity_cutoff){
			count++;
		}
	}
	return count;
}


void Sequences::add_usr_features(rulesTuplesList& feature_rules){
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


seqNames Sequences::get_names(){
  return m_sequence_names;
}


void Sequences::add_feature_indexes(FeaturesProfile& fprf){
  std::string nothing = "AA";
  for (auto &seq: m_sequences_aa){
    for (auto &res: seq){
        featureNamesList features = res.getFeatures();
        featuresList indexes;
        for (auto &feat: features){
          if (feat != nothing){
            indexes.push_back(fprf.findFeaturesIndex(feat));
          }
        }
        res.setFeatIndexes(indexes);
    }
  }
}
