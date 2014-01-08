//Profile class implementation
#include "Profile.h"
#include "substitutionMatrix.h"
#include "vecUtil.h"
#include "misc.h"
#include <string>
#include <iostream>
#include <vector>
Profile::Profile(std::vector< std::vector<double> > mat)
:	prfMatrix(mat){
}
Profile::Profile(){
}
Profile::Profile(Profile& that){
	prfMatrix = that.getMatrix();
}
Profile Profile::operator=(const Profile& that){
	Profile prfNew(that.getMatrix());
	return prfNew;
}
Profile::~Profile(){
	//cout << "DESTROYING PROFILE OBJECT\n";
}
//function getMatrix - returns profile matrix (double)
std::vector< std::vector<double> > Profile::getMatrix() const{
	return prfMatrix;
}
//calculate the alignment profile
void Profile::createProfile(std::vector<std::string>& alignment,const std::vector<bool>& sequenceIdentity){
	countOccurences(prfMatrix,alignment,sequenceIdentity);
	vecUtil::transposeVec(prfMatrix);
}
//calculate the alignment profile ENCODED SEQUENCES
void Profile::createProfile(std::vector< std::vector<std::string> >& alignment,const std::vector<bool>& sequenceIdentity, const std::vector<double>& sequenceIdentityValues, bool weightsModeOn){
	if (weightsModeOn) countOccurences(prfMatrix,alignment,sequenceIdentity);
	else countOccurences(prfMatrix,alignment,sequenceIdentity, sequenceIdentityValues);
	vecUtil::transposeVec(prfMatrix);
}
//builds a pseudo-profile from the profile itself and the substitution matrix with appropriate weights
void Profile::buildPseudoProfile(std::vector<std::string>& alignment, const std::vector<bool>& sequenceIdentity){
	createProfile(alignment,sequenceIdentity);
	std::vector< std::vector<double> > newProfile;
	for (int i = 0; i < prfMatrix[0].size(); i++){
		std::vector< std::vector<double> > columnsToAdd;
		for(int j = 0; j < prfMatrix.size(); j++){
			if (prfMatrix.at(j).at(i) != 0){
				std::vector<double> columnForJ = vecUtil::convertIntVectorToDoubleVector(substitutionMatrix::getColumn(j));
				vecUtil::multiplyVectorByAScalar(columnForJ, prfMatrix.at(j).at(i));
				columnsToAdd.push_back(columnForJ);
			}
		}
		newProfile.push_back(vecUtil::addUp(columnsToAdd));	//add up columns from substitution matrix for amino acids seen on ith position(times occurence/totalNrOfSeq))
	}
	vecUtil::transposeVec(newProfile);
	prfMatrix = newProfile;
}
//builds a pseudo-profile from the profile itself and the substitution matrix with appropriate weights for ENCODED SEQUENCES
void Profile::buildPseudoProfile(std::vector< std::vector<std::string> >& alignment, const std::vector<bool>& sequenceIdentity, const std::vector<double>& sequenceIdentityValues, bool weightsModeOn){
	createProfile(alignment,sequenceIdentity,sequenceIdentityValues,weightsModeOn);
	std::vector< std::vector<double> > newProfile;
	for (int i = 0; i < prfMatrix[0].size(); i++){
		std::vector< std::vector<double> > columnsToAdd;
		for(int j = 0; j < prfMatrix.size(); j++){
			if (prfMatrix.at(j).at(i) != 0){
				std::vector<double> columnForJ = vecUtil::convertIntVectorToDoubleVector(substitutionMatrix::getColumn(j));
				vecUtil::multiplyVectorByAScalar(columnForJ, prfMatrix.at(j).at(i));
				columnsToAdd.push_back(columnForJ);
			}
		}
		newProfile.push_back(vecUtil::addUp(columnsToAdd));	//add up columns from substitution matrix for amino acids seen on ith position(times occurence/totalNrOfSeq))
	}
	vecUtil::transposeVec(newProfile);
	prfMatrix = newProfile;
}
//function countOccurences returns matrix with occurences of each amino acid on each position normalized by the number of sequences
void Profile::countOccurences(std::vector< std::vector<double> >& result,std::vector<std::string>& alignment,const std::vector<bool>& sequenceIdentity){
	std::vector< std::vector<double> > tmpResult;
	int trueSequences = misc::countTrueValuesInVector(sequenceIdentity);
	for (int i = 0; i < alignment[0].size(); i++){
		std::vector<double> profileColumn(20,0);
		for (int j = 0; j < alignment.size(); j++){
			if (sequenceIdentity.at(j)){			//is it a sequence that I want to count in? (with identity > cutoff)
				char seqChar(alignment.at(j).at(i));
				if (seqChar != '-'){
					if (seqChar == 'B'){ 		//either D or N, so I'll add half a point to both
						profileColumn.at(2)+=0.5;
						profileColumn.at(3)+=0.5;
					}
					else if (seqChar == 'Z'){ 	//either D or N, so I'll add half a point to both
						profileColumn.at(6)+=0.5;
						profileColumn.at(7)+=0.5;
					}
					else if (seqChar == 'X'){
						for (int k = 0; k < profileColumn.size();k++){
							profileColumn.at(k)+=0.05;
						}
					}
					else{	
						int aAcidInt = substitutionMatrix::findAminoAcidsNo(seqChar);
						profileColumn.at(aAcidInt)++;				
					}
				}
			}
		}
		vecUtil::divideVectorByAScalar(profileColumn,trueSequences);
		tmpResult.push_back(profileColumn);
	}
	result = tmpResult;
}
//function countOccurences returns matrix with occurences of each amino acid on each position normalized by the number of sequences ENCODED SEQUENCES
void Profile::countOccurences(std::vector< std::vector<double> >& result,std::vector< std::vector<std::string> >& alignment,const std::vector<bool>& sequenceIdentity){
	std::vector< std::vector<double> > tmpResult;
	int trueSequences = misc::countTrueValuesInVector(sequenceIdentity);
	for (int i = 0; i < alignment[0].size(); i++){
		std::vector<double> profileColumn(20,0);
		for (int j = 0; j < alignment.size(); j++){
			if (sequenceIdentity.at(j)){			//is it a sequence that I want to count in? (with identity > cutoff)
				char seqChar(alignment.at(j).at(i)[0]);
				if (seqChar != '-'){
					if (seqChar == 'B'){ 		//either D or N, so I'll add half a point to both
						profileColumn.at(2)+=0.5;
						profileColumn.at(3)+=0.5;
					}
					else if (seqChar == 'Z'){ 	//either D or N, so I'll add half a point to both
						profileColumn.at(6)+=0.5;
						profileColumn.at(7)+=0.5;
					}
					else if (seqChar == 'X'){
						for (int k = 0; k < profileColumn.size();k++){
							profileColumn.at(k)+=0.05;
						}
					}
					else{	
						int aAcidInt = substitutionMatrix::findAminoAcidsNo(seqChar);
						profileColumn.at(aAcidInt)++;				
					}
				}
			}
		}
		vecUtil::divideVectorByAScalar(profileColumn,trueSequences);
		tmpResult.push_back(profileColumn);
	}
	result = tmpResult;
}
//function countOccurences returns matrix with occurences of each amino acid on each position normalized by the number of sequences ENCODED SEQUENCES with WEIGHTS
void Profile::countOccurences(std::vector< std::vector<double> >& result,std::vector< std::vector<std::string> >& alignment,const std::vector<bool>& sequenceIdentity, const std::vector<double>& sequenceIdentityValues){
	std::vector< std::vector<double> > tmpResult;
	double identitiesSum = vecUtil::sum(sequenceIdentityValues);
	for (int i = 0; i < alignment[0].size(); i++){
		std::vector<double> profileColumn(20,0);
		for (int j = 0; j < alignment.size(); j++){
			char seqChar(alignment.at(j).at(i)[0]);
			if (seqChar != '-'){
				if (seqChar == 'B'){ 		//either D or N, so I'll add half a point to both
					profileColumn.at(2)+=0.5*sequenceIdentityValues.at(j);
					profileColumn.at(3)+=0.5*sequenceIdentityValues.at(j);
				}
				else if (seqChar == 'Z'){ 	//either D or N, so I'll add half a point to both
					profileColumn.at(6)+=0.5*sequenceIdentityValues.at(j);
					profileColumn.at(7)+=0.5;
				}
				else if (seqChar == 'X'){
					for (int k = 0; k < profileColumn.size();k++){
						profileColumn.at(k)+=0.05*sequenceIdentityValues.at(j);
					}
				}
				else{	
					int aAcidInt = substitutionMatrix::findAminoAcidsNo(seqChar);
					profileColumn.at(aAcidInt)+=sequenceIdentityValues.at(j);				
				}
			}
		}
		vecUtil::divideVectorByAScalar(profileColumn,identitiesSum);
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
//function getElement - returns score for 'aAcid' amino acid on 'position' position
double Profile::getElement(int position, char aAcid){
	int result;
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
		int aAcidint = substitutionMatrix::findAminoAcidsNo(aAcid);
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
			std::cout << prfMatrix.at(i).at(j) << " ";
		}
		std::cout << "\n";
	}
}
//function printProfile - prints full profile
void Profile::printProfile(){
	for (int i = 0; i < prfMatrix.size(); i++){
		for (int j = 0; j < prfMatrix.at(0).size();j++){
			std::cout << prfMatrix.at(i).at(j) << " ";
		}
		std::cout << "\n";
	}
}
void Profile::setMatrix(std::vector<std::vector<double> > newMatrix){
	prfMatrix = newMatrix;
}
