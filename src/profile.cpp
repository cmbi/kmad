//Profile class implementation
#include "profile.h"
#include "residue.h"
#include "misc.h"
#include "substitution_matrix.h"
#include "vec_util.h"

#include<boost/range/numeric.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


Profile::Profile(ProfileMatrix mat)
:	m_prfMatrix(mat){
}


Profile::Profile(){
}


ProfileMatrix Profile::getMatrix() const{
	return m_prfMatrix;
}


void Profile::processProfile(SequenceList& alignment,
                             const IdentitiesList& sequenceIdentityValues){
	createProfile(alignment,sequenceIdentityValues);
	ProfileMatrix newProfile;
	for (unsigned int i = 0; i < m_prfMatrix[0].size(); i++){
		Matrix2D columnsToAdd;
		for(unsigned int j = 0; j < m_prfMatrix.size(); j++){
			if (m_prfMatrix[j][i] != 0){
        SbstMatColumn column_int; 
        substitutionMatrix::getColumn(j, column_int);
				ProfileMatrixColumn columnForJ = vecUtil::convertIntVectorToDoubleVector(column_int);
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


void Profile::createProfile(SequenceList& alignment, 
                            const IdentitiesList& sequenceIdentityValues){
	ProfileMatrix tmpResult;
	int noOfSequences = alignment.size();
	for (unsigned int i = 0; i < alignment[0].size(); i++){
		ProfileMatrixColumn profileColumn(20,0);
		int nonGaps = 0;
		for (unsigned int j = 0; j < alignment.size(); j++){
			char seqChar(alignment[j][i].getAA());
			if (seqChar != '-'){
				if (seqChar == 'B'){ 		//either D or N, so I'll add half a point to both
					profileColumn[2]+=0.5;
					profileColumn[3]+=0.5;
				}
				else if (seqChar == 'Z'){ 	//either D or N, so I'll add half a point to both
					profileColumn[6]+=0.5;
					profileColumn[7]+=0.5;
				}
				else if (seqChar == 'X'){
					for (unsigned int k = 0; k < profileColumn.size();k++){
						profileColumn[k]+=0.05;
					}
				}
				else{	
					int aAcidInt = substitutionMatrix::findAminoAcidsNo(seqChar);
					profileColumn[aAcidInt] += 1;				
				}
				nonGaps++;
			}
		}
		vecUtil::divideVectorByAScalar(profileColumn,noOfSequences);
		//vecUtil::divideVectorByAScalar(profileColumn,nonGaps);
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
