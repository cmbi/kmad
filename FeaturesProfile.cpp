#include "FeaturesProfile.h"
#include "vecUtil.h"
#include "misc.h"
#include <iostream>
#include <string>
#include <vector>
namespace {
	std::vector<std::string>listOfFeatures = {"phosphN", "domain0"};
}
FeaturesProfile::FeaturesProfile(std::vector< std::vector<double> > vec)
:	prfMatrix(vec){
}
FeaturesProfile::FeaturesProfile(){
}
void FeaturesProfile::createProfile(const std::vector< std::vector<std::string> >& alignment, const std::vector<bool>& sequenceIdentity, const std::vector<double>& sequenceIdentityValues, bool weightsModeOn){
	if (weightsModeOn) countOccurences(alignment, sequenceIdentityValues);
	else countOccurences(alignment,sequenceIdentity);
	vecUtil::transposeVec(prfMatrix);
}
//fucntion countOccurences - with weights mode on
void FeaturesProfile::countOccurences(const std::vector< std::vector<std::string> >& alignment, const std::vector<double>& sequenceIdentityValues){
	std::vector<std::vector<double> > tmpResult;
	double identitiesSum = vecUtil::sum(sequenceIdentityValues);
	char nothing = 'A';
	for (int i = 0; i < alignment[0].size();i++){
		std::vector<double> profileColumn(listOfFeatures.size(),0);
		for (int j = 0; j < alignment.size();j++){
				for(int k = 1; k < 4; k++){
					char alChar = alignment[j][i][k];
					if (alChar != nothing){
						int featIndex = findFeaturesIndex(name(alignment.at(j).at(i),k));
						if (featIndex != -1){
							if (k==3){
								profileColumn.at(featIndex)+=15*sequenceIdentityValues.at(j);
							}
							else if (k==2){
								profileColumn.at(featIndex)+=3*sequenceIdentityValues.at(j);
								profileColumn.at(1)-=3*sequenceIdentityValues.at(j);
							}
						}
					}
				}
		}
		vecUtil::divideVectorByAScalar(profileColumn,identitiesSum);
		tmpResult.push_back(profileColumn);
	} 
	prfMatrix = tmpResult;
}
void FeaturesProfile::countOccurences(const std::vector< std::vector<std::string> >& alignment, const std::vector<bool>& sequenceIdentity){
	std::vector<std::vector<double> > tmpResult;
	int trueSequences = misc::countTrueValuesInVector(sequenceIdentity);
	char nothing = 'A';
	for (int i = 0; i < alignment[0].size();i++){
		std::vector<double> profileColumn(listOfFeatures.size(),0);
		for (int j = 0; j < alignment.size();j++){
			if (sequenceIdentity.at(j)){
				for(int k = 1; k < 4; k++){
					char alChar = alignment[j][i][k];
					if (alChar != nothing){
						int featIndex = findFeaturesIndex(name(alignment.at(j).at(i),k));
						if (featIndex != -1){
							if (k==3){
								profileColumn.at(featIndex)+=15;
							}
							else if (k==2){
								profileColumn.at(featIndex)+=3;
								profileColumn.at(1)-=3;
							}
						}
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
	for (int i = 1; i < 4; i++){				//adds up scores from each feature (=each position from codon)
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
	}
	else if (featureType == 2){
		name="domain";
	}
	name.push_back(codon.at(featureType));
	return name;
}
//function expandListOfFeatures - expand it by domains and motifs found in the alignment
void FeaturesProfile::expandListOfFeatures(const std::vector< std::string >& sequence){
	std::string featureToAdd="";
	char nothing = 'A';
	for (int j = 0; j < sequence.size(); j++){
		char alChar = sequence.at(j).at(2);
		if (alChar != nothing){  										//means there is a domain on jth residue
			featureToAdd = "domain";
			featureToAdd += alChar;
			if (!vecUtil::contains(listOfFeatures,featureToAdd)){	//check whether this element is already in the list
				listOfFeatures.push_back(featureToAdd);
			}
		}
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
void FeaturesProfile::setMatrix(std::vector<std::vector<double> > newMatrix){
	prfMatrix = newMatrix;
}
