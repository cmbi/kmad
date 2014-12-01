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
:	m_prf_matrix(mat){
}


Profile::Profile(){
}


ProfileMatrix Profile::get_matrix() const{
	return m_prf_matrix;
}


void Profile::ProcessProfile(SequenceList& alignment){
	CreateProfile(alignment);
	ProfileMatrix new_profile;
	for (unsigned int i = 0; i < m_prf_matrix[0].size(); i++){
		Matrix2D columns_to_add;
		for(unsigned int j = 0; j < m_prf_matrix.size(); j++){
			if (m_prf_matrix[j][i] != 0){
        SbstMatColumn column_int; 
        substitutionMatrix::getColumn(j, column_int);
				ProfileMatrixColumn column_j = vecUtil::convertIntVectorToDoubleVector(column_int);
				vecUtil::multiplyVectorByAScalar(column_j, m_prf_matrix[j][i]);
				columns_to_add.push_back(column_j);
			}
		}
    //add up columns from substitution matrix for amino acids seen on ith 
    //position(times occurence/totalNrOfSeq))
		new_profile.push_back(vecUtil::addUp(columns_to_add));	
	}
	vecUtil::transposeVec(new_profile);
	m_prf_matrix = new_profile;
}


void Profile::CreateProfile(SequenceList& alignment){
	ProfileMatrix tmp_result;
	int no_of_sequences = alignment.size();
	for (unsigned int i = 0; i < alignment[0].size(); i++){
		ProfileMatrixColumn profile_column(20,0);
		int non_gaps = 0;
		for (unsigned int j = 0; j < alignment.size(); j++){
			char seq_char(alignment[j][i].getAA());
			if (seq_char != '-'){
				if (seq_char == 'B'){ 		//either D or N, so I'll add half a point to both
					profile_column[2]+=0.5;
					profile_column[3]+=0.5;
				}
				else if (seq_char == 'Z'){ 	//either D or N, so I'll add half a point to both
					profile_column[6]+=0.5;
					profile_column[7]+=0.5;
				}
				else if (seq_char == 'X'){
					for (unsigned int k = 0; k < profile_column.size();k++){
						profile_column[k]+=0.05;
					}
				}
				else{	
					int aacid_index = substitutionMatrix::findAminoAcidsIndex(seq_char);
					profile_column[aacid_index] += 1;				
				}
				non_gaps++;
			}
		}
		vecUtil::divideVectorByAScalar(profile_column,no_of_sequences);
		//vecUtil::divideVectorByAScalar(profileColumn,nonGaps);
		tmp_result.push_back(profile_column);
	}
	m_prf_matrix = tmp_result;
	vecUtil::transposeVec(m_prf_matrix);
}


double Profile::get_element(int position, char aacid){
	double result;
	if (aacid=='B'){ //take half the score for asparagine and half the score for aspartate
		result = 0.5*m_prf_matrix[2][position]+ 0.5*m_prf_matrix[3][position];
	}
	else if (aacid=='Z'){ // take half the score for glutamine and half the score for glutamate
		result = 0.5*m_prf_matrix[6][position]+ 0.5*m_prf_matrix[7][position];
	}
	else if (aacid=='X'){ // take average score from scores for all residues
		result = 0;
		for (unsigned int i = 0; i < m_prf_matrix.size(); i++){
			result += 0.05*m_prf_matrix[i][position];
		}
	}
	else { // it's not any of the {B,Z,X} -> single amino acid
		int aacid_index = substitutionMatrix::findAminoAcidsIndex(aacid);
		result = m_prf_matrix[aacid_index][position];
	}
	return result;
}


double Profile::get_element(int aacid_index, int position){
	return m_prf_matrix[aacid_index][position];
}
