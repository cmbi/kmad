#include "FeaturesProfile.h"
#include "vecUtil.h"
#include "txtProc.h"
#include "misc.h"
#include <iostream>
#include <string>
#include <vector>
namespace {
	std::vector<std::string> listOfFeatures = {"phosphN", "domain0", "motif0"};
}
//constructor, creates empty profile, takes domain and phosphorylation scores(dom, phosph) and motifs' ids and probabilities(m_ids, m_probs)
FeaturesProfile::FeaturesProfile(int dom, int phosph, std::vector<std::string> m_ids, std::vector<double> m_probs)
:	domainScore(dom),
	phosphScore(phosph),
	motifs_ids(m_ids),
	motifs_probs(m_probs)
	{
}
void FeaturesProfile::createProfile(const std::vector< std::vector<std::string> >& alignment, const std::vector<bool>& sequenceIdentity, const std::vector<double>& sequenceIdentityValues, bool weightsModeOn, int codon_length){
	if (weightsModeOn) countOccurences(alignment, sequenceIdentityValues, codon_length);
	else countOccurences(alignment,sequenceIdentity, codon_length);
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
					char alChar = alignment[j][i][k];
					if (alChar != nothing){
						int featIndex = findFeaturesIndex(name(alignment.at(j).at(i),k));
						if (featIndex != -1){
							if (k==3){
								profileColumn.at(featIndex)+=phosphScore * sequenceIdentityValues.at(j);
							}
							else if (k==2){
								profileColumn.at(featIndex) += domainScore * sequenceIdentityValues.at(j);
								profileColumn.at(1) -= domainScore * sequenceIdentityValues.at(j);
							}
						}
					}
					//has to be separate from domains and phosphorylations because it takes up 2 chars
					if (k == 4){
						std::string motifCode = "";
						motifCode.push_back(alignment[j][i][k]);
						motifCode.push_back(alignment[j][i][k+1]);
						int featIndex = findFeaturesIndex(name(alignment.at(j).at(i),k));
						double prob = motifs_prob(motifCode);
					}
				}
		}
		vecUtil::divideVectorByAScalar(profileColumn,identitiesSum);
		tmpResult.push_back(profileColumn);
	} 
	prfMatrix = tmpResult;
}
double FeaturesProfile::motifs_prob(std::string m_id){
	double prob;
	for (int i = 0; i < motifs_ids.size(); i++){
		if (motifs_ids[i] == m_id) {
			prob = motifs_probs[i];
			break;
		} 
	}
	return prob;
}
void FeaturesProfile::countOccurences(const std::vector< std::vector<std::string> >& alignment, const std::vector<bool>& sequenceIdentity, int codon_length){
	std::vector<std::vector<double> > tmpResult;
	int trueSequences = misc::countTrueValuesInVector(sequenceIdentity);
	char nothing = 'A';
	for (int i = 0; i < alignment[0].size();i++){
		std::vector<double> profileColumn(listOfFeatures.size(),0);
		for (int j = 0; j < alignment.size();j++){
			if (sequenceIdentity.at(j)){
				for(int k = 1; k < codon_length; k++){
					char alChar = alignment[j][i][k];
					if (alChar != nothing){
						int featIndex = findFeaturesIndex(name(alignment.at(j).at(i),k));
						if (featIndex != -1){
							if (k==3){
								profileColumn.at(featIndex) += phosphScore;
							}
							else if (k==2){
								profileColumn.at(featIndex) += domainScore;
								profileColumn.at(1)-= domainScore;
							}
						}
					}
					//has to be separate from domains and phosphorylations because it takes up 2 chars
					if (k == 4){
						std::string motifCode = "";
						motifCode.push_back(alignment[j][i][k]);
						motifCode.push_back(alignment[j][i][k+1]);
						int featIndex = findFeaturesIndex(name(alignment.at(j).at(i),k));
						double prob = motifs_prob(motifCode);


					}
				}
			}
		}
		vecUtil::divideVectorByAScalar(profileColumn,trueSequences);
		tmpResult.push_back(profileColumn);
	} 
	prfMatrix = tmpResult;
}
//function getScore - returns score for the entire codon on nth position
double FeaturesProfile::getScore(int position,std::string codon){
	double result = 0;
	for (int i = 1; i < codon.size(); i++){				//adds up scores from each feature (=each position from codon)
		result+=getElement(position,name(codon,i));	
	}
	return result;
}
//returns score for 1 feature on nth position
double FeaturesProfile::getElement(int position, std::string featName){
	int featuresIndex = findFeaturesIndex(featName);
	int result;
	if (featuresIndex == -1) result = 0;
	else result = prfMatrix.at(featuresIndex).at(position);	
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
	if (featureType == 3){
		name="phosph";
		name.push_back(codon.at(featureType));
	}
	else if (featureType == 2){
		name="domain";
		name.push_back(codon.at(featureType));
	}
	else if (featureType == 4){
		name="motif";
		name.push_back(codon.at(featureType));
		name.push_back(codon.at(featureType+1));
	}
	return name;
}
//function expandListOfFeatures - expand it by domains and motifs found in the alignment
void FeaturesProfile::expandListOfFeatures(const std::vector< std::string >& sequence, int codon_length){
	std::string featureToAdd="";
	char nothing = 'A';
	if (codon_length >= 3){
		for (int j = 0; j < sequence.size(); j++){
			char domainCode = sequence.at(j).at(2);
			if (domainCode != nothing){  										//means there is a domain on jth residue
				featureToAdd = "domain";
				featureToAdd += domainCode;
				if (!vecUtil::contains(listOfFeatures,featureToAdd)){	//check whether this element is already in the list
					listOfFeatures.push_back(featureToAdd);
				}
			}
		}
	}
	for (int i = 0; i < motifs_ids.size();i++){
		std::string motifCode = motifs_ids[i];
		featureToAdd = "motif";	
		featureToAdd += motifCode;
		listOfFeatures.push_back(featureToAdd);
	}
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
