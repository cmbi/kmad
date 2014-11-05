#include "FeaturesProfile.h"
#include "Residue.h"
#include "vecUtil.h"
#include "txtProc.h"
#include "misc.h"
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
namespace {
	// 0 - highest level of annotation, 3 - lowest, P - predicted
	std::vector<std::string> listOfFeatures = {"ptm_phosph0","ptm_phosph1","ptm_phosph2","ptm_phosph3","ptm_phosphP","ptm_acet0", "ptm_acet1","ptm_acet2","ptm_acet3","ptm_Nglyc0","ptm_Nglyc1","ptm_Nglyc2","ptm_Nglyc3","ptm_amid0","ptm_amid1","ptm_amid2","ptm_amid3","ptm_hydroxy0","ptm_hydroxy1","ptm_hydroxy2","ptm_hydroxy3","ptm_methyl0","ptm_methyl1","ptm_methyl2","ptm_methyl3","ptm_Oglyc0","ptm_Oglyc1","ptm_Oglyc2","ptm_Oglyc3","domain_0", "motif_0", "lcr"};
	std::vector<int> domain_indexes = {29};
	std::vector<int> motif_indexes = {30};
}
//constructor, creates empty profile, takes domain and phosphorylation scores(dom, phosph) and motifs' ids and probabilities(m_ids, m_probs), lcr - low complexity regions gap penlat modifier
FeaturesProfile::FeaturesProfile(int dom, int phosph, int motif, int lcr, 
                                 std::vector<std::string> m_ids, 
                                 std::vector<double> m_probs)
:	domainScore(dom),
	phosphScore(phosph),
	lcr_mod(lcr),
	motifScore(motif),
	motifs_ids(m_ids),
	motifs_probs(m_probs)
	{
}
// creates a features profile - counts occurences of each feature on each 
// position an normalize it by the number of non gaps (easily changeable to 
// number of sequences);(with or without weights)
void FeaturesProfile::createProfile(const std::vector< std::vector<Residue> >& alignment, 
                                    const std::vector<double>& sequenceIdentityValues, 
                                    bool weightsModeOn, int codon_length){
	std::vector<std::vector<double> > tmpResult;
	std::string nothing = "AA";
	double weight;
	double identitiesSum;
	int noOfSequences;
	if (weightsModeOn){
		identitiesSum = vecUtil::sum(sequenceIdentityValues);
	}
	else{
		noOfSequences = alignment.size();
	}
	for (unsigned int i = 0; i < alignment[0].size();i++){
		std::vector<double> profileColumn(listOfFeatures.size(),0);
		int nonGaps = 0;
		for (unsigned int j = 0; j < alignment.size();j++){
			if (alignment[j][i].getAA() != '-'){
				if (weightsModeOn){
					weight = sequenceIdentityValues[j];
				}
				else{
					weight = 1;
				}
				std::vector<std::string> features = alignment[j][i].getFeatures();	
				for (unsigned int k = 0; k < features.size(); k++){
					std::string featName = features[k];
					if (featName != nothing){
						int featIndex = findFeaturesIndex(featName);
						if (featIndex != -1){
							profileColumn[featIndex]+=modifier(featName) * weight; // seq. idenity * modifier (motif/domain/ptm...); w/o weights: weight == 1 
						}
						else{
							std::cout << "WARNING: unknown feature " << featName << std::endl;
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
		tmpResult.push_back(profileColumn);
	} 
	prfMatrix = tmpResult;
	vecUtil::transposeVec(prfMatrix);
}
//get motif's probability (by its id)
double FeaturesProfile::motifs_prob(std::string m_id){
	double prob = 0;
	std::string id_code = txtProc::split(m_id, '_')[1];
	for (unsigned int i = 0; i < motifs_ids.size(); i++){
		if (motifs_ids[i] == id_code) {
			prob = motifs_probs[i];
			break;
		} 
	}
	return prob;
}
//function getScore - returns score for the entire codon on nth position
void FeaturesProfile::getScore(int position, std::vector<std::string>& features, 
                               double& add_score, double& multiply_score, 
                               int sequence_no){
	std::string nothing = "AA";
	add_score = 0;
	multiply_score = 1;
	if (features[0] != nothing) score_PTMs(position,features[0], add_score, multiply_score);
	if (features[1] != nothing) score_domains(position,features[1], add_score, multiply_score);
	if (features[2] != nothing) score_motifs(position,features[2], add_score, multiply_score); 
	for (unsigned int i = 3; i < features.size(); i++){
		score_USR_features(sequence_no, position, features[i], add_score, multiply_score);	
	}
}
// return gap modifier, based on the low complexity regions on nth position. TO BE IMPLEMENTED
double FeaturesProfile::getGapMod(int position, std::vector<std::string> features){
	double result = 1;  // for now it doesn't change anything
	return result;
}
//returns score for a feature on nth position (by feature's name)
void FeaturesProfile::score_motifs(int position, std::string featName, double& add_score, double& multiply_score){
	int featuresIndex = findFeaturesIndex(featName);
	double result;
	if (featuresIndex == -1) result = 0;
	else result = prfMatrix[featuresIndex][position];	
	add_score += result;
}
//function findFeaturesIndex - takes features' name, e.g. "phosphN"
int FeaturesProfile::findFeaturesIndex(std::string featName){
	int featuresIndex = -1;
	std::string nothing = "AA";
	if (featName != nothing){
		for (unsigned int i = 0; i < listOfFeatures.size();i++){
			if (featName == listOfFeatures[i]){
				featuresIndex = i;
				break;
			}
		}
	}
	return featuresIndex;
}
//function name - converts codon to name, e.g. AAAN to phosphN
std::string FeaturesProfile::name(std::string codon, int featureType){
	std::string name;
	if (featureType == 4){
		switch(codon[featureType]){
			case 'N': name = "phosph0";
				break;
			case 'O': name = "phosph1"; 
				break;
			case 'P': name = "phosph2"; 
				break;
			case 'Q': name = "phosph3"; 
				break;
			case 'B': name = "acet0"; 
				break;
			case 'C': name = "acet1"; 
				break;
			case 'D': name = "acet2"; 
				break;
			case 'E': name = "acet3"; 
				break;
			case 'F': name = "Nglyc0"; 
				break;
			case 'G': name = "Nglyc1"; 
				break;
			case 'H': name = "Nglyc2"; 
				break;
			case 'I': name = "Nglyc3"; 
				break;
			case 'J': name = "amid0";
				break;
			case 'K': name = "amid1";
				break;
			case 'L': name = "amid2";
				break;
			case 'M': name = "amid3";
				break;
			case 'R': name = "hydroxy0";
				break;
			case 'S': name = "hydroxy1";
				break;
			case 'T': name = "hydroxy2";
				break;
			case 'U': name = "hydroxy3";
				break;
			case 'V': name = "methyl0";
				break;
			case 'W': name = "methyl1";
				break;
			case 'X': name = "methyl2";
				break;
			case 'Y': name = "methyl3";
				break;
			case 'Z': name = "Oglyc0";
				break;
			case 'a': name = "Oglyc1";
				break;
			case 'b': name = "Oglyc2";
				break;
			case 'c': name = "Oglyc3";
				break;
			case 'd': name = "phosphP";	//predicted phosphorylation
				break;
		}
		
	}
	else if (featureType == 2){
		name="domain";
		name.push_back(codon[featureType]);
		name.push_back(codon[featureType+1]);
	}
	else if (featureType == 5){
		name="motif";
		name.push_back(codon[featureType]);
		name.push_back(codon[featureType+1]);
	}
	else if (featureType == 1 && codon[featureType] == 'L'){
			name = "low_complexity_reg";
	}
	return name;
}
//function expandListOfFeatures - expand it by domains and motifs found in the alignment
void FeaturesProfile::expandListOfFeatures(const std::vector< std::vector<Residue> >& sequences){
	std::string nothing = "AA";  // code for no feature
	for(unsigned int i = 0; i < sequences.size();i++){	
		for (unsigned int j = 0; j < sequences[i].size(); j++){
			std::vector<std::string> features = sequences[i][j].getFeatures();
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
void FeaturesProfile::score_domains(int position, std::string dom_name, double& add_score, double& multiply_score){
	double result  = 0;
	std::string domain = "domain";
	for (unsigned int i = 0; i < domain_indexes.size(); i++){
		int dom_index = domain_indexes[i];
		if (listOfFeatures[dom_index] == dom_name){
			result += prfMatrix[dom_index][position];
		}
		else{
			result -= prfMatrix[dom_index][position];	// subtract score for every other domain
		}
	}
	add_score += result;
}
// get a score for a certain ptm encoded aligned to the 'position' position
void FeaturesProfile::score_PTMs(int position, std::string ptm_name, double& add_score, double& multiply_score){
	double result  = 0;
	std::string ptm_type = ptm_name;
	ptm_type.pop_back();	//pop back last character to get just the ptm type
	char ptm_level = ptm_name.back();	// level of annotation - last character of feature's name
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
	// now go through list of features to find in which rows in profile are the features that we're gonna score for
	for (unsigned int i = 0; i < listOfFeatures.size(); i++){
		std::string i_name = listOfFeatures[i];	
		std::string i_type = i_name;
		i_type.pop_back();  			// popping back last character, to get just the ptm type (without its level of annotation)
		if (i_type == ptm_type){
			char i_level = i_name.back();
			if (i_level == '0'){ result += prfMatrix[i][position];}
			else if (i_level == '1') result += prfMatrix[i][position] * 0.9;
			else if (i_level == '2') result += prfMatrix[i][position] * 0.8;
			else if (i_level == '3') result += prfMatrix[i][position] * 0.7;
			else if (i_level == 'P') result += prfMatrix[i][position] * 0.3;
		}
	}
	
	result = result * ptm_score;
	add_score += result;
}
void FeaturesProfile::score_USR_features(int sequence_no, int position, std::string feat_name, double& add_score, double& multiply_score){
	//first find tuple(s) with rules for this feature
	for (unsigned int i = 0; i < rules.size(); i++){
		if (std::get<0>(rules[i]) == feat_name && sequence_no == std::get<1>(rules[i])){
			double add_tmp = std::get<2>(rules[i]);
			double multiply_tmp = std::get<3>(rules[i]);
			std::vector<int> incr_features = std::get<6>(rules[i]);  //positions of increasing features
			for (unsigned int j = 0; j < incr_features.size(); j++){ //go through features that increase the score
				double prf_score = prfMatrix[incr_features[j]][position];
				if (prf_score != 0){
					add_score += add_tmp*prf_score;
					multiply_score *= 1 + (multiply_tmp - 1)*prf_score;
				}
			}
			//the same for decreasing features
			add_tmp = std::get<4>(rules[i]);
			multiply_tmp = std::get<5>(rules[i]);
			std::vector<int> decr_features = std::get<7>(rules[i]);
			for (unsigned int j = 0; j < decr_features.size(); j++){ //go through features that increase the score
				double prf_score = prfMatrix[decr_features[j]][position];
				if (prf_score != 0){
					add_score -= add_tmp*prf_score;
					multiply_score /= 1 + (multiply_tmp - 1)*prf_score;
				}
			}
		}
	}
}
void FeaturesProfile::printProfile(){
	for (unsigned int i = 0; i < prfMatrix.size(); i++){
		std::cout << listOfFeatures[i] << " ";
		for (unsigned int j = 0; j < prfMatrix[0].size();j++){
			std::cout << prfMatrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}
std::vector<std::vector<double> > FeaturesProfile::getMatrix(){
	return prfMatrix;
}
void FeaturesProfile::setMatrix(std::vector<std::vector<double> > newMatrix){
	prfMatrix = newMatrix;
}
//returns score modifier for a given feature
double FeaturesProfile::modifier(std::string featName){
	std::string feat_code = txtProc::split(featName,'_')[0];
	double modifier = 1;
	if (feat_code == "ptm"){
		modifier = phosphScore;
	}
	else if (feat_code == "domain"){
		modifier = domainScore;
	}
	else if (feat_code == "motif"){
		modifier = motifScore*motifs_prob(featName);
	}
	return modifier;
}


// sets rules for aligning features, unfolds last two elements in each tuple to vectors (vectors of features(strings))
void FeaturesProfile::setRules(std::vector<std::tuple<std::string,std::string,int,int,int,double,double,double,double,std::string,std::string> >& new_rules){
	for (unsigned int i = 0; i < new_rules.size(); i++){
		std::vector<int> incr_feat = txtProc::unfold(std::get<9>(new_rules[i]), listOfFeatures);
		std::vector<int> red_feat = txtProc::unfold(std::get<10>(new_rules[i]), listOfFeatures);
		// new_tuple <tag+name, incr_rule_1, incr_rule_2, red_rule_1, red_rule_2, incr_features_positions, red_features_positions>
		// incr_features_positions and red_features_positions are positions(indexes in the profile, not positions in sequence) of features that will be scored in the listOfFeatures vector
		std::tuple<std::string,int,double,double,double,double,std::vector<int>,std::vector<int>> new_tuple = std::make_tuple(std::string("USR_")+std::get<0>(new_rules[i])+std::string("_")+std::get<1>(new_rules[i]),
                                                                                                                          std::get<2>(new_rules[i]), 
                                                                                                                          std::get<5>(new_rules[i]),
                                                                                                                          std::get<6>(new_rules[i]),
                                                                                                                          std::get<7>(new_rules[i]), 
                                                                                                                          std::get<8>(new_rules[i]),
                                                                                                                          incr_feat, 
                                                                                                                          red_feat);
		rules.push_back(new_tuple);
	}
}


//adds user defined features fo listOfFeatures, called from txtProc::process_conf_file(....)
void FeaturesProfile::add_USR_features(std::vector<std::tuple<std::string, 
                                                              std::string,
                                                              int, int, int,
                                                              double, double,
                                                              double, double,
                                                              std::string,
                                                              std::string> >& new_rules){
	for (unsigned int i = 0; i < new_rules.size(); i++){
		std::string feature_i = std::string("USR_")+ std::get<0>(new_rules[i]) + std::string("_") + std::get<1>(new_rules[i]);
		if (!vecUtil::contains(listOfFeatures,feature_i)){
			listOfFeatures.push_back(feature_i);
		}
	}
}


void FeaturesProfile::printFeatures(){
	for (unsigned int i = 0; i < listOfFeatures.size(); i++){
		std::cout << listOfFeatures[i] << " ";
	}
	std::cout << std::endl;
}
