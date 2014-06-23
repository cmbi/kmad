#include "FeaturesProfile.h"
#include "vecUtil.h"
#include "txtProc.h"
#include "misc.h"
#include <iostream>
#include <string>
#include <vector>
namespace {
	std::vector<std::string> listOfFeatures = {"phosph0","phosph1","phosph2","phosph3","acet0", "acet1","acet2","acet3","Nglyc0","Nglyc1","Nglyc2","Nglyc3","amid0","amid1","amid2","amid3","hydroxy0","hydroxy1","hydroxy2","hydroxy3","methyl0","methyl1","methyl2","methyl3","Oglyc0","Oglyc1","Oglyc2","Oglyc3","domain0", "motif0", "low_complexity_reg"};
}
//constructor, creates empty profile, takes domain and phosphorylation scores(dom, phosph) and motifs' ids and probabilities(m_ids, m_probs), lcr - low complexity regions gap penlat modifier
FeaturesProfile::FeaturesProfile(int dom, int phosph, int motif, int lcr, std::vector<std::string> m_ids, std::vector<double> m_probs)
:	domainScore(dom),
	phosphScore(phosph),
	lcr_mod(lcr),
	motifScore(motif),
	motifs_ids(m_ids),
	motifs_probs(m_probs)
	{
}
void FeaturesProfile::createProfile(const std::vector< std::vector<std::string> >& alignment, const std::vector<double>& sequenceIdentityValues, bool weightsModeOn, int codon_length){
	if (weightsModeOn) countOccurences(alignment, sequenceIdentityValues, codon_length);
	else countOccurences(alignment, codon_length);
	vecUtil::transposeVec(prfMatrix);
}
//fucntion countOccurences - with weights mode on
void FeaturesProfile::countOccurences(const std::vector< std::vector<std::string> >& alignment, const std::vector<double>& sequenceIdentityValues, int codon_length){
	std::vector<std::vector<double> > tmpResult;
	double identitiesSum = vecUtil::sum(sequenceIdentityValues);
	char nothing = 'A';
	for (int i = 0; i < alignment[0].size();i++){
		std::vector<double> profileColumn(listOfFeatures.size(),0);
		for (int j = 0; j < alignment.size();j++){
				for(int k = 1; k < codon_length; k++){
					char alChar = alignment.at(j).at(i)[k];
					if (alChar != nothing && k != 3 && k != 6){  // 3rd and 6th position are for 2nd motif/domain character
						int featIndex = findFeaturesIndex(name(alignment.at(j).at(i),k));
						if (k==4 ){		//PTM (phosph)
							profileColumn.at(featIndex)+=phosphScore * sequenceIdentityValues.at(j);
						}
						else if (k==2){		//domain
							profileColumn.at(featIndex) += domainScore * sequenceIdentityValues.at(j);
						}
						else if (k == 5){	//motif
							std::string motifCode = "";
							motifCode.push_back(alignment[j][i][k]);
							motifCode.push_back(alignment[j][i][k+1]);
							double prob = motifs_prob(motifCode);
							profileColumn.at(featIndex) += motifScore * prob * sequenceIdentityValues.at(j);
						}
						else if (k==1){ 	//LCR
							profileColumn.at(featIndex) += 1;
						}
					}
				}
		}
		vecUtil::divideVectorByAScalar(profileColumn,identitiesSum);
		tmpResult.push_back(profileColumn);
	} 
	prfMatrix = tmpResult;
}
double FeaturesProfile::motifs_prob(std::string m_id){
	double prob = 0;
	for (int i = 0; i < motifs_ids.size(); i++){
		if (motifs_ids[i] == m_id) {
			prob = motifs_probs[i];
			break;
		} 
	}
	return prob;
}
void FeaturesProfile::countOccurences(const std::vector< std::vector<std::string> >& alignment, int codon_length){
	std::vector<std::vector<double> > tmpResult;
	int noOfSequences = alignment.size();
	std::string feature_code;
	char nothing = 'A';
	for (int i = 0; i < alignment[0].size();i++){
		std::vector<double> profileColumn(listOfFeatures.size(),0);
		for (int j = 0; j < alignment.size();j++){
			for(int k = 1; k < codon_length; k++){
				char alChar = alignment.at(j).at(i)[k];
				std::cout << alChar << " " << k << std::endl;
				if (alChar != nothing && k != 3 && k != 6){ // 3rd and 6th positions are for 2nd motif/domain character
					int featIndex = findFeaturesIndex(name(alignment.at(j).at(i),k));
					if (featIndex == -1 ) { std::cout << "that character ("<< alChar<<") shouldn't be on position "<< k<< std::endl;}
					if (k==3){
						profileColumn.at(featIndex) += phosphScore;
					}
					else if (k==2){
						profileColumn.at(featIndex) += domainScore;
					}
					else if (k == 4){
						feature_code = "";
						feature_code.push_back(alignment[j][i][k]);
						feature_code.push_back(alignment[j][i][k+1]);
						double prob = motifs_prob(feature_code);
						profileColumn.at(featIndex) += motifScore * prob;
					}
					else if (k==1){
						profileColumn.at(featIndex) += 1;
					}
				}
			}
		}
		vecUtil::divideVectorByAScalar(profileColumn,noOfSequences);
		tmpResult.push_back(profileColumn);
	} 
	prfMatrix = tmpResult;
}
//function getScore - returns score for the entire codon on nth position
double FeaturesProfile::getScore(int position,std::string codon){
	double result = 0;
	char nothing = 'A';
	for (int i = 2; i < codon.size()-1; i++){	//adds up scores from each feature (=each position from codon, except for the first(0) one - amino acid and the second one - lcr, which influences only the gap penalty)
		if (nothing != codon[i] && i != 2 && i != 3 && i != 6 ){	//positions 3 and 6 a second character of motif/domain code, domains scored separately (because of subtraction for 'other' domains)
			result+=getElement(position,name(codon,i));
			if (i == 2){
				result += score_domains(position,name(codon,i));
			}
		}
	}
	return result;
}
double FeaturesProfile::getGapMod(int position, std::string codon){
	double result = 1;  // for now it doesn't change anything
	/*
	if (codon[1]=='L'){ // if it's a low complexity region in the sequence
		result = 0.5*lcr_mod;
	}
	else{
		result  = 0.5;
	}
	double lcr_profile = getElement(position,"low_complexity_reg");//lcr frequency at this position 
	result += (1 - (1 - lcr_mod)*lcr_profile)*0.5; // score from profile ranges from 0.5*lcr_modifier to 0.5, depending on the frequency in the profile 
	*/
	return result;
}
//returns score for 1 feature on nth position
double FeaturesProfile::getElement(int position, std::string featName){
	int featuresIndex = findFeaturesIndex(featName);
	int result;
	if (featuresIndex == -1) result = 0;
	else result = prfMatrix[featuresIndex][position];	
	return result;
}
double FeaturesProfile::getElement(int position, int featuresIndex){
	return prfMatrix.at(featuresIndex).at(position);
}
//function findFeaturesIndex - takes features' name, e.g. "phosphN"
int FeaturesProfile::findFeaturesIndex(std::string featName){
	int featuresIndex = -1;
	for (int i = 0; i < listOfFeatures.size();i++){
		if (featName == listOfFeatures[i]){
			featuresIndex = i;
			break;
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
		}
		
	}
	else if (featureType == 2){
		name="domain";
		name.push_back(codon.at(featureType));
		name.push_back(codon.at(featureType+1));
	}
	else if (featureType == 5){
		name="motif";
		name.push_back(codon.at(featureType));
		name.push_back(codon.at(featureType+1));
	}
	else if (featureType == 1 && codon[featureType] == 'L'){
			name = "low_complexity_reg";
	}
	return name;
}
//function expandListOfFeatures - expand it by domains and motifs found in the alignment
void FeaturesProfile::expandListOfFeatures(const std::vector< std::string >& sequence, int codon_length){
	std::string featureToAdd="";
	std::string nothing = "AA";
	std::string feature_code;
	if (codon_length >= 3){
		for (int j = 0; j < sequence.size(); j++){
			feature_code = "";
			feature_code.push_back(sequence.at(j).at(2));
			feature_code.push_back(sequence.at(j).at(3));
			if (feature_code != nothing){  										//means there is a domain on jth residue
				featureToAdd = "domain";
				featureToAdd += feature_code;
				if (!vecUtil::contains(listOfFeatures,featureToAdd)){	//check whether this domain is already in the list of features
					listOfFeatures.push_back(featureToAdd);
				}
			}
		}
	}
	for (int i = 0; i < motifs_ids.size();i++){    // i can expand motifs' list differently than domains, because motifs ids are already stored in "motifs_ids"
		std::string motifCode = motifs_ids[i];
		featureToAdd = "motif";	
		featureToAdd += motifCode;
		listOfFeatures.push_back(featureToAdd);
	}
}
double FeaturesProfile::score_domains(int position, std::string dom_name){
	double result  = 0;
	std::string domain = "domain";
	for (int i = 0; i < listOfFeatures.size(); i++){
		std::string i_name = listOfFeatures[i];	
		i_name.pop_back();
		i_name.pop_back();
		if (i_name == domain){
			if (listOfFeatures[i] == dom_name){
				result += prfMatrix[i][position];
			}
			else{
				result -= prfMatrix[i][position];
			}
		}
	}
	return result;
}
void FeaturesProfile::printProfile(){
	for (int i = 0; i < prfMatrix.size(); i++){
		for (int j = 0; j < prfMatrix.at(0).size();j++){
			std::cout << prfMatrix.at(i).at(j) << " ";
		}
		std::cout << std::endl;
	}
}
std::vector<std::vector<double> > FeaturesProfile::getMatrix(){
	return prfMatrix;
}
std::vector<std::string> FeaturesProfile::getMotifsIDs(){
	return motifs_ids;
}
std::vector<double> FeaturesProfile::getMotifsProbs(){
	return motifs_probs;
}
double FeaturesProfile::getDomainScore(){
	return domainScore;
}
double FeaturesProfile::getPhosphScore(){
	return phosphScore;
}
void FeaturesProfile::setMatrix(std::vector<std::vector<double> > newMatrix){
	prfMatrix = newMatrix;
}
