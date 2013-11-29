#include "FeaturesProfile.h"
#include "vecUtil.h"
#include <iostream>
#include <string>
#include <vector>
namespace {
	std::vector<std::string>listOfFeatures = {"phosphN"};
}
FeaturesProfile::FeaturesProfile(std::vector< std::vector<double> > vec)
:	prfMatrix(vec){
}
FeaturesProfile::FeaturesProfile(){
}
void FeaturesProfile::createProfile(const std::vector< std::vector<std::string> >& alignment, const std::vector<bool>& sequenceIdentity){
	expandListOfFeatures(alignment,sequenceIdentity);
	countOccurences(alignment,sequenceIdentity);
	vecUtil::transposeVec(prfMatrix);
}
void FeaturesProfile::countOccurences(const std::vector< std::vector<std::string> >& alignment, const std::vector<bool>& sequenceIdentity){
	std::vector<std::vector<double> > tmpResult;
	for (int i = 0; i < alignment[0].size();i++){
		std::vector<double> profileColumn(listOfFeatures.size(),0);
		for (int j = 0; j < alignment.size();j++){
			if (sequenceIdentity.at(j)){
				for(int k = 1; k < 4; k++){
					profileColumn[findFeaturesIndex(name(alignment.at(j).at(i),k))];
				}
			}
		}
	} 
}
//function getScore - returns score for the entire codon on nth position
double FeaturesProfile::getScore(std::string codon, int position){
	double result = 0;
	for (int i = 1; i < 4; i++){
		result+=getElement(position,name(codon,i));	
	}
	return result;
}
//returns score for 1 feature on nth position
double FeaturesProfile::getElement(int position, std::string featName){
	int featuresIndex = findFeaturesIndex(featName);
	return prfMatrix.at(featuresIndex).at(position);	
}
double FeaturesProfile::getElement(int position, int featuresIndex){
	return prfMatrix.at(featuresIndex).at(position);
}
int FeaturesProfile::findFeaturesIndex(std::string featName){
	int featuresIndex;
	for (int i = 0; i < listOfFeatures.size();i++){
		if (featName == listOfFeatures[i]){
			featuresIndex = i;
			break;
		}
	}
	return featuresIndex;
}
std::string FeaturesProfile::name(std::string codon, int featureType){
	std::string name;
	if (featureType == 2){
		name="phosph";
	}
	else if (featureType == 1){
		name="domain";
	}
	name.push_back(codon.at(featureType));
	return name;
}
//function expandListOfFeatures - expand it by domains and motifs found in alignment
void FeaturesProfile::expandListOfFeatures(const std::vector< std::vector<std::string> >& alignment, const std::vector<bool>& sequenceIdentity){
/*
	string featureToAdd;
	for (int i = 0; i < sequenceIdentity.size(); i++){
		if (sequenceIdentity[i]){
			for (int j = 0; j < alignment.at(i).size(); j++){
				for (int k = 1; k < 4; k++){
					if (alignment.at(i).at(k) != "A"){
						if (k==3){		//PTMs are coded by 3rd char
							if (alignment.at(i).at(k)=="N"){
								featureToAdd = "phosphN";
							}
						}
					}
				}
			}
		}
	}
*/
}

