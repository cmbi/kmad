#include <stdio>
#include <string>
#include <vector>
#include "FeaturesProfile.h"
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
void FeaturesProfile::countOccurences(const vector< vector<string>& alignment >, const vector<bool>& sequenceIdentity){
	for (int i = 0; i < alignment.size();i++){
		for (int j = 0; j < alignment.at(i).size();j++){
			
		}
	} 
}
double FeaturesProfile::getElement(int position, string features){
	char phosph = features[3];
	int featuresIndex = findFeaturesIndex(phosph, "phosph");
	return prfMatrix.at(featuresIndex).at(position);	
}
double FeaturesProfile::getElement(int position, int featuresIndex){
	return prfMatrix.at(featuresIndex).at(position);
}
int FeaturesProfile::findFeaturesIndex(char feature,string featureType){
	feat = featureType;
	feat.push_back(feature);
	int featuresIndex;
	for (int i = 0; i < listOfFeatures.size();i++){
		if (feat == listOfFeatures[i]){
			featuresIndex = i;
			break;
		}
	}
	return featuresIndex;
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

