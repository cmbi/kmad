//Profile class implementation
#include <string>
#include <iostream>
#include <vector>
#include "Profile.h"
#include "SubstitutionMatrix.h"
#include "UsefulStuff.h"
using namespace std;
extern UsefulStuff util;
Profile::Profile(vector< vector<double> > vec){
	prfMatrix = vec;
}
Profile::Profile(){
}
Profile::Profile(Profile& that){
	prfMatrix = that.getMatrix();
}
/*Profile& operator=(Profile& that){
	prfMatrix = that.getMatrix();
}*/
Profile::~Profile(){
	//cout << "DESTROYING PROFILE OBJECT\n";
}
//calculate the alignment profile
void Profile::createProfile(vector<string>& alignment,const vector<bool>& sequenceIdentity){
	countOccurences(prfMatrix,alignment,sequenceIdentity);
	util.transposeVec(prfMatrix);
}
//builds a pseudo-profile from the profile itself and the substitution matrix with appropriate weights
void Profile::buildPseudoProfile(vector<string>& alignment, const vector<bool>& sequenceIdentity, SubstitutionMatrix& sbst){
	createProfile(alignment,sequenceIdentity);
	vector< vector<double> > newProfile;
	for (int i = 0; i < prfMatrix[0].size(); i++){
		vector< vector<double> > columnsToAdd;
		for(int j = 0; j < prfMatrix.size(); j++){
			if (prfMatrix.at(j).at(i) != 0){
				vector<double> columnForJ = util.convertIntVectorToDoubleVector(sbst.getColumn(j));
				util.multiplyVectorByAScalar(columnForJ, prfMatrix.at(j).at(i));
				columnsToAdd.push_back(columnForJ);
			}
		}
		newProfile.push_back(util.addUp(columnsToAdd));	//add up columns from substitution matrix for amino acids seen on ith position(times occurence/totalNrOfSeq))
	}
	util.transposeVec(newProfile);
	prfMatrix = newProfile;
}
//function countOccurences returns matrix with occurences of each amino acid on each position normalized by the number of sequences
//vector< vector<double> > 
void Profile::countOccurences(vector< vector<double> >& result,vector<string>& alignment,const vector<bool>& sequenceIdentity){
	vector< vector<double> > tmpResult;
	int trueSequences = util.countTrueValuesInVector(sequenceIdentity);
	for (int i = 0; i < alignment[0].size(); i++){
		vector<double> profileColumn(20,0);
		for (int j = 0; j < alignment.size(); j++){
			if (sequenceIdentity.at(j)){			//is it a sequence that I want to count in? (with identity > 80%)
				char seqChar(alignment.at(j).at(i));
				if (seqChar != '-'){
					if (seqChar == 'B'){ 		//either D or N, so I'll add half a point to 
						profileColumn.at(2)+=0.5;
						profileColumn.at(3)+=0.5;
					}
					else if (seqChar == 'Z'){ 	//either D or N, so I'll add half a point to 
						profileColumn.at(6)+=0.5;
						profileColumn.at(7)+=0.5;
					}
					else if (seqChar == 'X'){
						for (int k = 0; k < profileColumn.size();k++){
							profileColumn.at(k)+=0.05;
						}
					}
					else{	
						int aAcidInt = util.findAminoAcidsNo(seqChar);
						profileColumn.at(aAcidInt)++;				
					}
				}
			}
		}
		util.divideVectorByAScalar(profileColumn,trueSequences);
		tmpResult.push_back(profileColumn);
	}
	result = tmpResult;
}
//function countNonGaps - counts how many characters in alignment nth ('column' integer) column are not gaps 
double Profile::countNonGaps(int column){
	double sum = 0;
	for (int i =0; i < 20; i++){
		sum += double(prfMatrix.at(i).at(column));
	}
	return sum;
}

//function getMatrix - returns profile matrix (double)
vector< vector<double> > Profile::getMatrix(){
	return prfMatrix;
}
//function getElement - returns score for 'aAcid' amino acid on 'position' position
double Profile::getElement(int position, char aAcid){
	double result;
	if (aAcid=='B'){
		result = 0.5*prfMatrix.at(2).at(position)+ 0.5*prfMatrix.at(3).at(position);
	}
	else if (aAcid=='Z'){
		result = 0.5*prfMatrix.at(6).at(position)+ 0.5*prfMatrix.at(7).at(position);
	}
	else if (aAcid=='X'){
		result = 0;
		for (int i = 0; i < prfMatrix.size(); i++){
			result += 0.05*prfMatrix.at(i).at(position);
		}
	}
	else {	
		int aAcidint = util.findAminoAcidsNo(aAcid);
		result = prfMatrix.at(aAcidint).at(position);
	}
	return result;
}
//function getElement - returns score for nth amino acid on mth position
double Profile::getElement(int aAcidInt, int position){
	return prfMatrix.at(aAcidInt).at(position);
}
//function printProfile(int,int) - prints only columns from boundStart to boundEnd
void Profile::printProfile(int boundStart, int boundEnd){
	for (int i = boundStart; i < boundEnd; i++){
		for (int j = 0; j < prfMatrix.at(0).size();j++){
			cout << prfMatrix.at(i).at(j) << " ";
		}
		cout << endl;
	}
}
//function printProfile - prints full profile
void Profile::printProfile(){
	for (int i = 0; i < prfMatrix.size(); i++){
		for (int j = 0; j < prfMatrix.at(0).size();j++){
			cout << prfMatrix.at(i).at(j) << " ";
		}
		cout << endl;
	}
}
