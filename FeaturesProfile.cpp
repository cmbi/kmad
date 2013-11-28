#include <iostream>
#include <string>
#include <vector>
#include "FeaturesProfile.h"
#include "UsefulStuff.h"
using namespace std;
extern UsefulStuff util;
FeaturesProfile::FeaturesProfile(vector< vector<double> > vec){
	prfMatrix = vec;
	listOfFeatures.push_back("phosphN");
}
FeaturesProfile::FeaturesProfile(){
}
void FeaturesProfile::createProfile(const vector< vector<string> >& alignment, const vector<bool>& sequenceIdentity){
	expandListOfFeatures(alignment,sequenceIdentity);
	countOccurences(alignment,sequenceIdentity);
	util.transposeVec(prfMatrix);
}
void FeaturesProfile::countOccurences(const vector< vector<string> >& alignment, const vector<bool>& sequenceIdentity){
	vector<vector<double> > tmpResult;
	for (int i = 0; i < alignment[0].size();i++){
		vector<double> profileColumn(listOfFeatures.size(),0);
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
double FeaturesProfile::getScore(string codon, int position){
	double result = 0;
	for (int i = 1; i < 4; i++){
		result+=getElement(position,name(codon,i));	
	}
	return result;
}
//returns score for 1 feature on nth position
double FeaturesProfile::getElement(int position, string featName){
	int featuresIndex = findFeaturesIndex(featName);
	return prfMatrix.at(featuresIndex).at(position);	
}
double FeaturesProfile::getElement(int position, int featuresIndex){
	return prfMatrix.at(featuresIndex).at(position);
}
int FeaturesProfile::findFeaturesIndex(string featName){
	int featuresIndex;
	for (int i = 0; i < listOfFeatures.size();i++){
		if (featName == listOfFeatures[i]){
			featuresIndex = i;
			break;
		}
	}
	return featuresIndex;
}
string FeaturesProfile::name(string codon, int featureType){
	string name;
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
void FeaturesProfile::expandListOfFeatures(const vector< vector<string> >& alignment, const vector<bool>& sequenceIdentity){
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

