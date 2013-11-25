//SubtitutionMatrix class implementation
#include <iostream>
#include <vector>
#include "SubstitutionMatrix.h"
#include "UsefulStuff.h"
using namespace std;
extern UsefulStuff util;
//constructor
SubstitutionMatrix::SubstitutionMatrix(){
	static const int simScoresArr[20][20] = { { 4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0},{-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3},{-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3},{-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3},{ 0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},{-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2},{-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2},{ 0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3},{-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3},{-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3},{-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1},{-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2},{-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1},{-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1},{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2},{ 1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2},{ 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0},{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3},{-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1},{ 0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4} };
	static const char alphabetArr[] = { 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V' };
	for (int i=0; i < 20;i++){
		alphabet.push_back(alphabetArr[i]);
	}
	vector<int> newEntry(20);
	for (int i=0; i < 20;i++){
		for (int j = 0; j < 20;j++){
			newEntry.at(j) = simScoresArr[i][j];
		}
		simScores.push_back(newEntry);	
	}
}
//function getElement - returns score for substitution of char1 amino acid by char2
int SubstitutionMatrix::getElement(char char1,char char2){
	int index1=0;
	int index2=0;
	for (int i = 0; i<20;i++){
		if(alphabet.at(i)==char1){
			index1 = i;
		}
		if(alphabet.at(i)==char2){
			index2 = i;
		}
	}
	return simScores.at(index1).at(index2);
}
//convert to profile format - creates a matrix 20 x sequence length, where nth column is a column from sbstMatrix for the amino acid on position n in the sequence
vector< vector<double> > SubstitutionMatrix::convertToProfileFormat(string sequence){
	vector< vector<double> > result(sequence.size());
	for (int i = 0; i < result.size(); i++){
		int aAcidInt = util.findAminoAcidsNo(sequence[i]);
		result.at(i)=util.convertIntVectorToDoubleVector(simScores.at(aAcidInt));//adds a column to the result(converted from int to double)
	}
	util.transposeVec(result);
	return result;
}
//function getElement - retruns score for substitution of ith aa by jth aa
int SubstitutionMatrix::getElement(int i,int j){
	return simScores.at(i).at(j);
}
//function printSbstMatrix - prints substitution matrix
void SubstitutionMatrix::printSbstMatrix(){
	for (int i=0; i < simScores.size();i++){
		for (int j = 0; j < simScores.at(i).size();j++){
			cout << simScores.at(i).at(j);
			cout << " ";
		}
		cout << "\n";
	}
}
//function getLetter - returns nth('lNr' integer) aa
char SubstitutionMatrix::getLetter(int lNr){
	return alphabet.at(lNr);
}
vector<int> SubstitutionMatrix::getColumn(int columnNo){
	return simScores.at(columnNo);
}
