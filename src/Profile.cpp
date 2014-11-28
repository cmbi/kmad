//Profile class implementation
#include "Profile.h"
#include "Residue.h"
#include "misc.h"
#include "substitutionMatrix.h"
#include "vecUtil.h"

#include<boost/range/numeric.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


Profile::Profile(profile_matrix mat)
:	m_prfMatrix(mat){
}


Profile::Profile(){
}


//function getMatrix - returns profile matrix (double)
profile_matrix Profile::getMatrix() const{
	return m_prfMatrix;
}


//builds a pseudo-profile from the profile itself and the substitution matrix with appropriate weights
void Profile::buildPseudoProfile(sequenceList& alignment,
                                 const identitiesList& sequenceIdentityValues, 
                                 bool weightsModeOn){
	createProfile(alignment,sequenceIdentityValues,weightsModeOn);
	profile_matrix newProfile;
	for (unsigned int i = 0; i < m_prfMatrix[0].size(); i++){
		matrix2d columnsToAdd;
		for(unsigned int j = 0; j < m_prfMatrix.size(); j++){
			if (m_prfMatrix[j][i] != 0){
        sbstMatColumn column_int; 
        substitutionMatrix::getColumn(j, column_int);
				profileMatrixColumn columnForJ = vecUtil::convertIntVectorToDoubleVector(column_int);
				vecUtil::multiplyVectorByAScalar(columnForJ, m_prfMatrix[j][i]);
				columnsToAdd.push_back(columnForJ);
			}
		}
    //add up columns from substitution matrix for amino acids seen on ith 
    //position(times occurence/totalNrOfSeq))
		newProfile.push_back(vecUtil::addUp(columnsToAdd));	
	}
	vecUtil::transposeVec(newProfile);
	m_prfMatrix = newProfile;
}


void Profile::createProfile(sequenceList& alignment, 
                            const identitiesList& sequenceIdentityValues,
                            bool weightsModeOn){
	profile_matrix tmpResult;
	double identitiesSum;
	int noOfSequences;
	if (weightsModeOn){
    identitiesSum = boost::accumulate(sequenceIdentityValues, 0);
	}
	else{
		noOfSequences = alignment.size();
	}
	for (unsigned int i = 0; i < alignment[0].size(); i++){
		profileMatrixColumn profileColumn(20,0);
		int nonGaps = 0;
		for (unsigned int j = 0; j < alignment.size(); j++){
			
			double weight;
			if (weightsModeOn){
				weight = sequenceIdentityValues[j];
			}
			else{
				weight = 1;
			}
			char seqChar(alignment[j][i].getAA());
			if (seqChar != '-'){
				if (seqChar == 'B'){ 		//either D or N, so I'll add half a point to both
					profileColumn[2]+=0.5*weight;
					profileColumn[3]+=0.5*weight;
				}
				else if (seqChar == 'Z'){ 	//either D or N, so I'll add half a point to both
					profileColumn[6]+=0.5*weight;
					profileColumn[7]+=0.5*weight;
				}
				else if (seqChar == 'X'){
					for (unsigned int k = 0; k < profileColumn.size();k++){
						profileColumn[k]+=0.05*weight;
					}
				}
				else{	
					int aAcidInt = substitutionMatrix::findAminoAcidsNo(seqChar);
					profileColumn[aAcidInt]+=weight;				
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
	m_prfMatrix = tmpResult;
	vecUtil::transposeVec(m_prfMatrix);
}


//function getElement - returns score for 'aAcid' amino acid on 'position' position
double Profile::getElement(int position, char aAcid){
	double result;
	if (aAcid=='B'){ //take half the score for asparagine and half the score for aspartate
		result = 0.5*m_prfMatrix[2][position]+ 0.5*m_prfMatrix[3][position];
	}
	else if (aAcid=='Z'){ // take half the score for glutamine and half the score for glutamate
		result = 0.5*m_prfMatrix[6][position]+ 0.5*m_prfMatrix[7][position];
	}
	else if (aAcid=='X'){ // take average score from scores for all residues
		result = 0;
		for (unsigned int i = 0; i < m_prfMatrix.size(); i++){
			result += 0.05*m_prfMatrix[i][position];
		}
	}
	else { // it's not any of the {B,Z,X} -> single amino acid
		int aAcidint = substitutionMatrix::findAminoAcidsNo(aAcid);
		result = m_prfMatrix[aAcidint][position];
	}
	return result;
}


//function getElement - returns score for nth amino acid on mth position
double Profile::getElement(int aAcidInt, int position){
	return m_prfMatrix[aAcidInt][position];
}
