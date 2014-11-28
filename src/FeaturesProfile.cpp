#include "FeaturesProfile.h"
#include "Residue.h"
#include "vecUtil.h"
#include "txtProc.h"
#include "misc.h"

#include<boost/range/numeric.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
namespace {
	// 0 - highest level of annotation, 3 - lowest, P - predicted
	featureNamesList listOfFeatures = {"ptm_phosph0", "ptm_phosph1",
                                     "ptm_phosph2", "ptm_phosph3",
                                     "ptm_phosphP", "ptm_acet0",
                                     "ptm_acet1", "ptm_acet2", 
                                     "ptm_acet3", "ptm_Nglyc0",
                                     "ptm_Nglyc1", "ptm_Nglyc2",
                                     "ptm_Nglyc3", "ptm_amid0",
                                     "ptm_amid1", "ptm_amid2",
                                     "ptm_amid3", "ptm_hydroxy0",
                                     "ptm_hydroxy1", "ptm_hydroxy2",
                                     "ptm_hydroxy3", "ptm_methyl0",
                                     "ptm_methyl1", "ptm_methyl2",
                                     "ptm_methyl3", "ptm_Oglyc0",
                                     "ptm_Oglyc1", "ptm_Oglyc2",
                                     "ptm_Oglyc3", "domain_0",
                                     "motif_0", "lcr"};
	indexList domain_indexes = {29};
	indexList motif_indexes = {30};
  std::string nothing = "AA";
	std::string domain = "domain";
}


//constructor, creates empty profile, takes domain and phosphorylation 
//scores(dom, phosph) and motifs' ids and probabilities(m_ids, m_probs), 
//lcr - low complexity regions gap penlat modifier
FeaturesProfile::FeaturesProfile(int dom, int phosph, int motif, int lcr, 
                                 ids_list m_ids, 
                                 probs_list m_probs)
:	m_domainScore(dom),
	m_phosphScore(phosph),
	m_motifScore(motif),
	m_motifs_ids(m_ids),
	m_motifs_probs(m_probs)
	{
}


// creates a features profile - counts occurences of each feature on each 
// position an normalize it by the number of non gaps (easily changeable to 
// number of sequences);(with or without weights)
void FeaturesProfile::createProfile(const sequenceList& alignment, 
                                    const identitiesList& sequenceIdentityValues, 
                                    bool weightsModeOn, int codon_length){
  countOccurences(alignment, sequenceIdentityValues, 
                  weightsModeOn, codon_length);
  processProfile();
}
void FeaturesProfile::processProfile(){
  m_prfMatrix.clear();
  for (unsigned int i = 0; i < m_occurences_matrix.size(); i++){
    std::string feat_name = listOfFeatures[i];
    profileMatrixRow feature_row(m_occurences_matrix[i].size());
    for (unsigned int j = 0; j < m_occurences_matrix[i].size(); j++){
      if (feat_name.substr(0,3) == "ptm"){
        feature_row[j] = score_PTMs(j, feat_name);
      }
      else if (feat_name.substr(0,6) == "domain"){
        feature_row[j] = score_domains(j, feat_name);
      }
      else if (feat_name.substr(0,5) == "motif"){
        feature_row[j] = score_motifs(j, feat_name);
      }
      else if (feat_name.substr(0,3) == "USR"){
        feature_row[j] = score_USR_features(j, feat_name);
      }
    }
    m_prfMatrix.push_back(feature_row);
  }
}
void FeaturesProfile::countOccurences(const sequenceList& alignment, 
                                    const identitiesList& sequenceIdentityValues, 
                                    bool weightsModeOn, int codon_length){
  m_occurences_matrix.clear();
	std::string nothing = "AA";
	double weight;
	double identitiesSum;
	int noOfSequences;
	if (weightsModeOn){
    identitiesSum = boost::accumulate(sequenceIdentityValues, 0);
	}
	else{
		noOfSequences = alignment.size();
	}
	for (unsigned int i = 0; i < alignment[0].size();i++){
		profileMatrixColumn profileColumn(listOfFeatures.size(),0);
		int nonGaps = 0;
		for (unsigned int j = 0; j < alignment.size();j++){
			if (alignment[j][i].getAA() != '-'){
				if (weightsModeOn){
					weight = sequenceIdentityValues[j];
				}
				else{
					weight = 1;
				}
				featureNamesList features = alignment[j][i].getFeatures();	
				for (unsigned int k = 0; k < features.size(); k++){
					std::string feat_name = features[k];
					if (feat_name != nothing){
						int featIndex = findFeaturesIndex(feat_name);
						if (featIndex != -1){
              // seq. idenity * modifier (motif/domain/ptm...); 
              // w/o weights: weight == 1 
							profileColumn[featIndex] += get_modifier(feat_name) * weight; 
						}
					}
				}
			nonGaps++;
			}
		}
		if (weightsModeOn){
			vecUtil::divideVectorByAScalar(profileColumn,identitiesSum);
		}
		else{
			vecUtil::divideVectorByAScalar(profileColumn,noOfSequences);
			//vecUtil::divideVectorByAScalar(profileColumn,nonGaps);
    }
		m_occurences_matrix.push_back(profileColumn);
	} 
	vecUtil::transposeVec(m_occurences_matrix);
}
//get motif's probability (by its id)
double FeaturesProfile::get_motifs_prob(std::string& m_id){
	double prob = 0;
	std::string id_code = txtProc::split(m_id, '_')[1];
	for (unsigned int i = 0; i < m_motifs_ids.size(); i++){
		if (m_motifs_ids[i] == id_code) {
			prob = m_motifs_probs[i];
			break;
		} 
	}
	return prob;
}


void FeaturesProfile::getScore(unsigned int position, featuresList& features, 
                               double& add_score){
  for (unsigned int i = 0; i < features.size(); i++){
       add_score += m_prfMatrix[features[i]][position];
  }
}


//returns score for a feature on nth position (by feature's name)
double FeaturesProfile::score_motifs(unsigned int& position, std::string& feat_name){
	int featuresIndex = findFeaturesIndex(feat_name);
	double result;
	if (featuresIndex == -1) result = 0;
	else result = m_occurences_matrix[featuresIndex][position];	
  return result;
}


//function findFeaturesIndex - takes features' name, e.g. "phosphN"
int FeaturesProfile::findFeaturesIndex(std::string& feat_name){
	int featuresIndex = -1;
	for (unsigned int i = 0; i < listOfFeatures.size();i++){
		if (feat_name == listOfFeatures[i]){
			featuresIndex = i;
			break;
		}
	}
	return featuresIndex;
}


//function expandListOfFeatures - expand it by domains and motifs found in the alignment
void FeaturesProfile::expandListOfFeatures(const sequenceList& sequences){
	for(unsigned int i = 0; i < sequences.size();i++){	
		for (unsigned int j = 0; j < sequences[i].size(); j++){
			featureNamesList features = sequences[i][j].getFeatures();
			for (unsigned int k = 0; k < features.size();k++){
				std::string feature_k = features[k];
				if (!vecUtil::contains(listOfFeatures,feature_k) && feature_k != nothing){	//check whether this domain is already in the list of features
					listOfFeatures.push_back(feature_k);
					if (feature_k.substr(0,6) == "domain"){
						domain_indexes.push_back(listOfFeatures.size()-1); //to later look for domains only in these positions (not necessary but saves time)
					}
					else if (feature_k.substr(0,5) == "motif"){
						motif_indexes.push_back(listOfFeatures.size()-1); // same as with domain_indexes
					}
				}
			}
		}
	}
}


//get score for domain "dom_name" on certain position
double FeaturesProfile::score_domains(unsigned int& position, std::string& dom_name){
	double result  = 0;
	for (unsigned int i = 0; i < domain_indexes.size(); i++){
		int dom_index = domain_indexes[i];
		if (listOfFeatures[dom_index] == dom_name){
			result += m_occurences_matrix[dom_index][position];
		}
		else{
			result -= m_occurences_matrix[dom_index][position];	// subtract score for every other domain
		}
	}
  return result;
}


// get a score for a certain ptm encoded aligned to the 'position' position
double FeaturesProfile::score_PTMs(unsigned int& position, std::string& ptm_name){
	double result  = 0;
	std::string ptm_type = ptm_name;
  //pop back last character to get just the ptm type
	ptm_type.pop_back();	
  // level of annotation - last character of feature's name
	char ptm_level = ptm_name.back();	
	double ptm_score;
	// first set ptm_score based on annotation level of the query ptm
	if (ptm_level == '0'){ ptm_score = 1.0;}
	else if (ptm_level == '1'){ ptm_score = 0.9;}
	else if (ptm_level == '2'){ ptm_score = 0.8;}
	else if (ptm_level == '3'){ ptm_score = 0.7;}
	else if (ptm_level == 'P'){ ptm_score = 0.3;}
	else {
		std::cout << "something's wrong with annotation level on position "<< position << " ptmname: " << ptm_name<< std::endl;
		std::exit(0);
	}
	// now go through list of features to find in which rows in profile 
  // are the features that we're gonna score for
	for (unsigned int i = 0; i < listOfFeatures.size(); i++){
		std::string i_name = listOfFeatures[i];	
		std::string i_type = i_name;
    // popping back last character, to get just the ptm type 
    // (without its level of annotation)
		i_type.pop_back();  			
		if (i_type == ptm_type){
			char i_level = i_name.back();
			if (i_level == '0'){ result += m_occurences_matrix[i][position];}
			else if (i_level == '1') result += m_occurences_matrix[i][position] * 0.9;
			else if (i_level == '2') result += m_occurences_matrix[i][position] * 0.8;
			else if (i_level == '3') result += m_occurences_matrix[i][position] * 0.7;
			else if (i_level == 'P') result += m_occurences_matrix[i][position] * 0.3;
		}
	}
	result = result * ptm_score;
  return result;
}


double FeaturesProfile::score_USR_features(unsigned int& position, 
                                           std::string& feat_name){
	//first find tuple(s) with rules for this feature
  double result = 0;
  for (auto &rule: m_rules){
		if (std::get<0>(rule) == feat_name){
			double add_tmp = std::get<2>(rule);
      //positions of increasing features
			featuresList incr_features = std::get<6>(rule);  
      //go through features that increase the score
      for (auto &feat: incr_features){
				double prf_score = m_occurences_matrix[feat][position];
				if (prf_score != 0){
					 result += add_tmp*prf_score;
				}
			}
			//the same for decreasing features
			add_tmp = std::get<4>(rule);
			featuresList decr_features = std::get<7>(rule);
      //go through features that increase the score
			//for (unsigned int j = 0; j < decr_features.size(); j++){ 
      for (auto &feat: decr_features){
				double prf_score = m_occurences_matrix[feat][position];
				if (prf_score != 0){
					result -= add_tmp*prf_score;
				}
			}
		}
	}
  return result;
}
 

profile_matrix FeaturesProfile::getMatrix(){
	return m_prfMatrix;
}


//returns score modifier for a given feature
double FeaturesProfile::get_modifier(std::string& feat_name){
	std::string feat_code = txtProc::split(feat_name,'_')[0];
	double modifier = 1;
	if (feat_code == "ptm"){
		modifier = m_phosphScore;
	}
	else if (feat_code == "domain"){
		modifier = m_domainScore;
	}
	else if (feat_code == "motif"){
		modifier = m_motifScore*get_motifs_prob(feat_name);
	}
	return modifier;
}


// sets rules for aligning features, unfolds last two elements in each tuple 
// to vectors (vectors of features(strings))
void FeaturesProfile::setRules(rulesTuplesList& new_rules){
  for (auto &rule: new_rules){
		featuresList incr_feat = txtProc::unfold(std::get<9>(rule), listOfFeatures);
		featuresList red_feat = txtProc::unfold(std::get<10>(rule), listOfFeatures);
		// new_tuple <tag+name, incr_rule_1, incr_rule_2, red_rule_1, red_rule_2, 
    // incr_features_positions, red_features_positions>
		// incr_features_positions and red_features_positions are positions(indexes
    // in the profile, not positions in sequence) of features that will be 
    // scored in the listOfFeatures vector
		processedRules new_tuple = std::make_tuple(std::string("USR_") 
                                                + std::get<0>(rule) 
                                                + std::string("_")
                                                + std::get<1>(rule),
                                               std::get<2>(rule), 
                                               std::get<5>(rule),
                                               std::get<6>(rule),
                                               std::get<7>(rule), 
                                               std::get<8>(rule),
                                               incr_feat, 
                                               red_feat);
		m_rules.push_back(new_tuple);
	}
}


//adds user defined features fo listOfFeatures, called from 
//txtProc::process_conf_file(....)
void FeaturesProfile::add_USR_features(rulesTuplesList& new_rules){
  for (auto &rule: new_rules){
		std::string feature_i = std::string("USR_") + std::get<0>(rule) \
                            + std::string("_") + std::get<1>(rule);
		if (!vecUtil::contains(listOfFeatures,feature_i)){
			listOfFeatures.push_back(feature_i);
		}
	}
}
