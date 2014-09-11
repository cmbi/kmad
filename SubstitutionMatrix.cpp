//SubtitutionMatrix class implementation
#include "substitutionMatrix.h"
#include "misc.h"
#include "vecUtil.h"
#include <iostream>
#include <vector>
//constructor
namespace {
	static const std::vector<char>alphabet = { 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V' };
	//BLOSUM62
	static const std::vector< std::vector<int> >simScores = {{4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},{-1,-3, -3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3},{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1},{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}};
	//DISORDER
	//static const std::vector< std::vector<int> > simScores = {{3, -2, -1, -1, -1, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -5, -2, 0}, {-2, 5, 0, -2, -1, 1, -1, -2, 0, -2, -2, 2, -1, -3, -2, -1, -1, 0, -2, -2}, {-1, 0, 4, 1, -1, 1, 0, 0, 2, -3, -3, 0, -2, -2, -1, 1, 0, -3, -1, -3}, {-1, -2, 1, 4, -3, 0, 2, -1, -1, -4, -4, -1, -4, -4, -2, 0, -1, -4, -4, -4}, {-1, -1, -1, -3, 10, -3, -4, -3, -1, 0, -1, -3, 0, -1, -2, 0, 1, -5, 0, 1}, {-1, 1, 1, 0, -3, 5, 0, -2, 1, -2, -2, 0, -1, -2, -1, 0, 0, -1, 0, -2}, {-1, -1, 0, 2, -4, 0, 4, -2, -1, -3, -3, 0, -3, -4, -1, -1, -1, -4, -3, -2}, {0, -2, 0, -1, -3, -2, -2, 5, -1, -5, -4, -2, -4, -4, -1, 0, -2, -4, -3, -4}, {-2, 0, 2, -1, -1, 1, -1, -1, 8, -2, -2, -1, -2, 0, -2, -1, 0, -2, 2, -2}, {-1, -2, -3, -4, 0, -2, -3, -5, -2, 4, 2, -2, 1, 1, -2, -2, -1, -2, 0, 3}, {-1, -2, -3, -4, -1, -2, -3, -4, -2, 2, 4, -2, 2, 1, -1, -2, -2, -2, 0, 1}, {-1, 2, 0, -1, -3, 0, 0, -2, -1, -2, -2, 4, -2, -3, -1, -1, 0, -3, -2, -2}, {-1, -1, -2, -4, 0, -1, -3, -4, -2, 1, 2, -2, 7, 1, -2, -2, -1, -1, -1, 1}, {-2, -3, -2, -4, -1, -2, -4, -4, 0, 1, 1, -3, 1, 7, -3, -2, -2, -1, 4, 0}, {-1, -2, -1, -2, -2, -1, -1, -1, -2, -2, -1, -1, -2, -3, 6, 0, -1, -1, -3, -1}, {1, -1, 1, 0, 0, 0, -1, 0, -1, -2, -2, -1, -2, -2, 0, 3, 1, -3, -2, -2}, {0, -1, 0, -1, 1, 0, -1, -2, 0, -1, -2, 0, -1, -2, -1, 1, 4, -5, -1, 0}, {-5, 0, -3, -4, -5, -1, -4, -4, -2, -2, -2, -3, -1, -1, -1, -3, -5, 13, 3, -4}, {-2, -2, -1, -4, 0, 0, -3, -3, 2, 0, 0, -2, -1, 4, -3, -2, -1, 3, 8, -1}, {0, -2, -3, -4, 1, -2, -2, -4, -2, 3, 1, -2, 1, 0, -1, -2, 0, -4, -1, 4}};
	//dummy 0 5
	//static const std::vector< std::vector<int> > simScores = {{5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5}};
	//DUNKER extended alphabet	
	//static const std::vector< std::vector<int> > simScores = {{9, 1, 2, 2, 1, 3, 2, 2, 1, 1, 1, 2, 0, 0, 3, 3, 4, -1, 0, 3}, {1, 10, 1, 0, 1, 3, 1, 2, 3, 0, 0, 5, -1, -2, 0, 1, 1, 2, -1, 0}, {2, 1, 11, 4, 1, 3, 2, 2, 3, 0, -1, 3, -1, -1, 1, 4, 3, -2, 1, 0}, {2, 0, 4, 10, -2, 2, 5, 2, 2, -1, -2, 1, -2, -2, 1, 2, 2, -3, -1, 0}, {1, 1, 1, -2, 17, 0, -3, 1, 0, 0, 0, -1, -2, 1, -1, 2, 1, 3, 2, 1}, {3, 3, 3, 2, 0, 11, 4, 0, 4, 0, 1, 3, 0, -1, 2, 2, 2, 1, 0, 1}, {2, 1, 2, 5, -3, 4, 9, 1, 1, -1, -1, 3, -1, -2, 0, 1, 1, -2, -2, 0}, {2, 2, 2, 2, 1, 0, 1, 10, 0, -2, -2, 0, -2, -2, 0, 2, 1, 0, -2, 0}, {1, 3, 3, 2, 0, 4, 1, 0, 13, 0, 0, 1, -1, 2, 1, 1, 1, 1, 4, 0}, {1, 0, 0, -1, 0, 0, -1, -2, 0, 12, 5, 0, 4, 4, 0, 0, 2, 1, 2, 7}, {1, 0, -1, -2, 0, 1, -1, -2, 0, 5, 10, -1, 4, 5, 1, 0, 1, 2, 2, 4}, {2, 5, 3, 1, -1, 3, 3, 0, 1, 0, -1, 10, -1, -2, 0, 1, 2, -2, -1, 0}, {0, -1, -1, -2, -2, 0, -1, -2, -1, 4, 4, -1, 13, 2, -2, -1, 1, 0, 0, 3}, {0, -2, -1, -2, 1, -1, -2, -2, 2, 4, 5, -2, 2, 13, -1, 0, 0, 6, 8, 3}, {3, 0, 1, 1, -1, 2, 0, 0, 1, 0, 1, 0, -2, -1, 11, 2, 2, -1, -2, 1}, {3, 1, 4, 2, 2, 2, 1, 2, 1, 0, 0, 1, -1, 0, 2, 9, 4, -1, 0, 1}, {4, 1, 3, 2, 1, 2, 1, 1, 1, 2, 1, 2, 1, 0, 2, 4, 10, -2, 0, 3}, {-1, 2, -2, -3, 3, 1, -2, 0, 1, 1, 2, -2, 0, 6, -1, -1, -2, 18, 6, 0}, {0, -1, 1, -1, 2, 0, -2, -2, 4, 2, 2, -1, 0, 8, -2, 0, 0, 6, 14, 1}, {3, 0, 0, 0, 1, 1, 0, 0, 0, 7, 4, 0, 3, 3, 1, 1, 3, 0, 1, 11}};
	//dummy binary
	//static const std::vector<std::vector<int> > simScores = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}};
}
//function getElement - returns score for substitution of char1 amino acid by char2
int substitutionMatrix::getElement(char char1,char char2){
	int index1=0;
	int index2=0;
	for (int i = 0; i<20;i++){
		if(alphabet[i] == char1){
			index1 = i;
		}
		if(alphabet[i] == char2){
			index2 = i;
		}
	}
	return simScores[index1][index2];
}
/*
//convert to profile format - creates a matrix 20 x sequence length, where nth column is a column from sbstMatrix for the amino acid on position n in the sequence
std::vector< std::vector<double> > substitutionMatrix::convertToProfileFormat(std::string sequence){
	std::vector< std::vector<double> > result(sequence.size());
	for (int i = 0; i < result.size(); i++){
		int aAcidInt = findAminoAcidsNo(sequence[i]);
		result.at(i)=vecUtil::convertIntVectorToDoubleVector(simScores.at(aAcidInt));//adds a column to the result(converted from int to double)
	}
	vecUtil::transposeVec(result);
	return result;
}
*/
//convert to profile format - creates a matrix 20 x sequence length, where nth column is a column from sbstMatrix for the amino acid on position n in the sequence ENCODED SEQUENCES
std::vector< std::vector<double> > substitutionMatrix::convertToProfileFormat(std::vector<std::string> sequence){
	std::vector< std::vector<double> > result(sequence.size());
	std::vector<std::vector<int> > newSbstRow;
	for (int i = 0; i < result.size(); i++){
		if (sequence[i][0] == 'B'){
			newSbstRow.clear();
			newSbstRow.push_back(simScores[2]);
			newSbstRow.push_back(simScores[3]);
			result[i] = vecUtil::average(newSbstRow);
		}
		else if (sequence[i][0] == 'Z'){
			newSbstRow.clear();
			newSbstRow.push_back(simScores[6]);
			newSbstRow.push_back(simScores[7]);
			result[i] = vecUtil::average(newSbstRow);

		}
		else if (sequence[i][0] == 'X'){
			result[i] = vecUtil::average(simScores);
		}
		else{
			int aAcidInt = findAminoAcidsNo(sequence[i][0]);
			//std:: cout << aAcidInt << " " << sequence[i][0] << std::endl;
			result[i] = vecUtil::convertIntVectorToDoubleVector(simScores[aAcidInt]);//adds a column to the result(converted from int to double)
		}
	}
	vecUtil::transposeVec(result);
	return result;
}
//function getElement - retruns score for substitution of ith aa by jth aa
int substitutionMatrix::getElement(int i,int j){
	return simScores[i][j];
}
//function printSbstMatrix - prints substitution matrix
void substitutionMatrix::printSbstMatrix(){
	for (int i=0; i < simScores.size();i++){
		for (int j = 0; j < simScores[i].size();j++){
			std::cout << simScores[i][j];
			std::cout << " ";
		}
		std::cout << "\n";
	}
}
//function getLetter - returns nth('lNr' integer) aa
char substitutionMatrix::getLetter(int lNr){
	return alphabet[lNr];
}
std::vector<int> substitutionMatrix::getColumn(int columnNo){
	return simScores[columnNo];
}
//function findAminoAcidsNo - finds index of the given char aa amino acid
int substitutionMatrix::findAminoAcidsNo(char aa){
	int aAcidint = -1;
	for (int i = 0; i < alphabet.size();i++){
		if (aa == alphabet[i]){
			aAcidint = i;
			break;
		}
	}
	return aAcidint;
}
