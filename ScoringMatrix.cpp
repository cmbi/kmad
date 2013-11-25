//ScoringMatrix class implementation
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "ScoringMatrix.h"
#include "SubstitutionMatrix.h"
#include "Profile.h"
#include "UsefulStuff.h"
using namespace std;
extern UsefulStuff util;
/*constructor
	arguments:
		s1size - length of the 1st sequence
		s2size - length of the 2nd sequence
		pen - gap opening penalty
*/
//constructor
ScoringMatrix::ScoringMatrix(int s1size,int s2size, int pen){
	iLength = s1size;
	jLength = s2size;	
	gapOpening = double(pen);
	gapExtension = double(-1);
	gapOpeningHorizontal = gapOpening;
	gapExtensionHorizontal = gapExtension;
	vector<double> row(jLength+1,0);
	matrixV.assign(iLength+1, row);
	matrixG.assign(iLength+1, row);
	matrixH.assign(iLength+1, row);
}
//destructor
ScoringMatrix::~ScoringMatrix(){
//	cout <<  "DESTROYING SCORING MATRIX\n";
}
ScoringMatrix::ScoringMatrix(ScoringMatrix& that){
	iLength = that.iLength;		
	jLength = that.jLength;
	gapOpening = that.gapOpening;
	gapExtension = that.gapExtension;
	gapOpeningHorizontal = gapOpening;
	gapExtensionHorizontal = gapExtension;
	matrixV = that.matrixV;
	matrixG = that.matrixG;
	matrixH = that.matrixH;
}
//function calculateScoresProfile - calculates scoring matrix for sequences s1 and s2 using profile prf instead of a substitution matrix
void ScoringMatrix::calculateScores(string s2, Profile& prf, int debug){
	string s1(prf.getMatrix().at(0).size()+1,'A');	//makes a pseudosequence of the length of the profile+1 (poly-A)
	s2 = string("0").append(s2);
	string s1String,s2String;
	for (int i = 1; i < matrixV.size(); i++){
		matrixV.at(i).at(0) = -10000000; //infinity
		matrixH.at(i).at(0) = -10000000;
		matrixG.at(i).at(0) = 0;	//makes gaps at the beginning of the vertical sequence (s1) free
		//matrixG[i][0] = gapOpening+(i-1)*gapExtension;
	}
	for (int i = 1; i < matrixV[0].size(); i++){
		matrixV.at(0).at(i) = -10000000;
		matrixH.at(0).at(i) = 0;	//makes gaps at the end of the horizontal sequence (s2) free
		//matrixH[0][i] = gapOpening+(i-1)*gapExtension;
		matrixG.at(0).at(i) = -10000000;
	}
	double score1,score2,score3;
	for (int i = 1; i < matrixV.size();i++){
		for (int j = 1; j < matrixV.at(i).size(); j++){
			///V
			score1 = matrixV.at(i-1).at(j-1) + prf.getElement(i-1,s2.at(j));
			score2 = matrixG.at(i-1).at(j-1) + prf.getElement(i-1,s2.at(j));
			score3 = matrixH.at(i-1).at(j-1) + prf.getElement(i-1,s2.at(j));
			matrixV[i][j] = util.maxValue(score1,score2,score3);
			///G
			score1 = matrixV.at(i-1).at(j) + gapOpening;
			score2 = matrixG.at(i-1).at(j) + gapExtension;
			matrixG[i][j] = (score1 > score2) ? score1 : score2;
			///H
			score1 = matrixV.at(i).at(j-1) + gapOpeningHorizontal;
			score2 = matrixH.at(i).at(j-1) + gapExtensionHorizontal;
			matrixH[i][j] = (score1 > score2) ? score1 : score2;
		}
	}
}
//function findBestScore - returns alignment score with positions in the scoring matrix: [score, i, j] (must be either in the last row or in the last column of the scoring matrix)
vector<int> ScoringMatrix::findBestScore(){
	int maxI = matrixV.size()-1;
	int maxJ = matrixV.at(0).size()-1;
	int n = maxI;
	int m = maxJ;
	double maxIval = matrixV.at(maxI).at(maxJ);
	double maxJval = matrixV.at(maxI).at(maxJ);
	double max = matrixV.at(maxI).at(maxJ);
	for (int i = 0; i < n ; i++){		//finds max score in the last row
		if (matrixV.at(i).at(m) > max){
			maxIval = matrixV.at(i).at(m);
			max = matrixV.at(i).at(m);
			maxI = i;
		} 
	}
	for (int i = 0; i < m; i++){		//finds max score in the last column	
		if (matrixV.at(n).at(i)>max){
			maxJval = matrixV.at(n).at(i);
			max = matrixV.at(n).at(i);
			maxJ = i;
		}
	}
	vector<int> resArr;
	if (maxIval > maxJval){			//max score is in the last row
		resArr.push_back(maxI);
		resArr.push_back(m);
	}
	else{					//max score is in the last column
		resArr.push_back(n);		// max score position (i)
		resArr.push_back(maxJ);		//max score position (j)
	}
	return resArr;
}
//function getVec - returns scoring matrix
vector< vector<double> > ScoringMatrix::getVec(){
	return matrixV;
}
//function nwAlignment - performs a sequence vs profile(/pseudoprofile) needleman wunsch alignment
//vector<string> 
void ScoringMatrix::nwAlignment(vector<string> *result,string s2, Profile& prf,bool verbose){
	string s1(prf.getMatrix()[0].size()+1,'A');
	s2 = string("0").append(s2);
	string newS1 = "";
	string newS2 = "";
	vector<string> ali;
	string newChar1="";	
	string newChar2="";	
	int i = s1.length()-1;
	int j = s2.length()-1;
	string currentMatrix = "V";
	int iteratorM = 10;
	//if bestScore isn't in the lower right corner, then add gaps to newS1 or newS2
	if (findBestScore().at(0) != matrixV.size()-1 || findBestScore().at(1) != matrixV.at(0).size()-1){
		i = findBestScore().at(0);
		j = findBestScore().at(1);
		for (int k = s1.length()-1; k > i; k--){
			newChar1 = s1.at(k);
			newChar2 = "-";
			newS1 = newChar1.append(newS1);
			newS2 = newChar2.append(newS2);
		}
		for (int k = s2.length()-1; k > j; k--){
			newChar1 = "-";
			newChar2 = s2.at(k);
			newS2 = newChar2.append(newS2);
			newS1 = newChar1.append(newS1);
		}
	}
	//trace back the matrix
	while (i > 0 || j > 0){
		if (i > 0 && j > 0 && currentMatrix == "V"){	//match/mismatch
			newChar1 = s1.at(i);
			newChar2 = s2.at(j);
			if (matrixV.at(i).at(j) != matrixV.at(i-1).at(j-1) + prf.getElement(i-1,s2.at(j))){
				if( i > 0 && j > 0 && matrixV.at(i).at(j) == matrixG.at(i-1).at(j-1)+prf.getElement(i-1,s2.at(j))){
					currentMatrix = "G";	
				}
				else if(i > 0 && j > 0 && matrixV.at(i).at(j) == matrixH.at(i-1).at(j-1)+prf.getElement(i-1,s2.at(j))){
					currentMatrix = "H";
				}
			}
			i--;
			j--;
		}
		else if (i > 0 && currentMatrix == "G"){	//gap in seq2
			newChar1 = s1.at(i);
			newChar2 = "-";
			if (matrixG.at(i).at(j) == matrixV.at(i-1).at(j) + gapOpening)
				currentMatrix = "V";
			i--;
		}
		else if (j > 0 && currentMatrix == "H"){	//gap in profile
			newChar1 = "-";
			newChar2 = s2.at(j);
			if (matrixH.at(i).at(j) == matrixV.at(i).at(j-1) + gapOpeningHorizontal){
				currentMatrix = "V";
			}
			j--;
		}
		newS1 = newChar1.append(newS1);
		newS2 = newChar2.append(newS2);
	}
	if (verbose){
		cout << newS1 << endl << newS2 << endl << endl;
	}
	ali.push_back(newS1);
	ali.push_back(newS2);
	*result = ali;
}
