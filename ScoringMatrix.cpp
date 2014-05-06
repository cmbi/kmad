//ScoringMatrix class implementation
//
//responsible for pairwise alignments
//creates scoring matrix, backtraces alignment path, etc.
#include "ScoringMatrix.h"
#include "FeaturesProfile.h"
#include "substitutionMatrix.h"
#include "Profile.h"
#include "findVal.h"
#include "misc.h"
#include "vecUtil.h"
#include "txtProc.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
/*constructor
	arguments:
		s1size - length of the 1st sequence
		s2size - length of the 2nd sequence
		pen - gap opening penalty
*/
//constructor
ScoringMatrix::ScoringMatrix(int s1size,int s2size, double pen, double extensionPenalty)
:	iLength(s1size),
	jLength(s2size),	
	gapOpening(double(pen)),
	gapExtension(double(extensionPenalty)),
	gapOpeningHorizontal(gapOpening),
	gapExtensionHorizontal(gapExtension)
{
	std::vector<double> row(jLength+1,0);
	matrixV.assign(iLength+1, row);
	matrixG.assign(iLength+1, row);
	matrixH.assign(iLength+1, row);
}
//function calculateScoresProfile - calculates scoring matrix for sequences s1 and s2 using profile prf instead of a substitution matrix ENCODED SEQUENCES
void ScoringMatrix::calculateScores(std::vector<std::string> s2, Profile& prf, FeaturesProfile& featPrf,int debug, int codon_length){
	std::vector<std::string> s1 = misc::pseudoSequence(prf.getMatrix()[0].size()+1); //creating polyA pseudoSequence representing the profile, to know later where are the gaps in the profile
	s2 = vecUtil::push_front(s2,txtProc::gapCode(codon_length));
	std::string s1String,s2String;
	for (int i = 1; i < matrixV.size(); i++){
		matrixV[i][0] = -10000000; //infinity
		matrixH[i][0] = -10000000;
		matrixG[i][0] = 0;	//makes gaps at the beginning of the vertical sequence (s1) free
		//matrixG[i][0] = gapOpening+(i-1)*gapExtension; // uncomment to penalize gaps at the beginnig of s1 seq
	}
	for (int i = 1; i < matrixV[0].size(); i++){
		matrixV[0][i] = -10000000;
		matrixH[0][i] = 0;	//makes gaps at the end of the horizontal sequence (s2) free
		//matrixH[0][i] = gapOpening+(i-1)*gapExtension; //uncomment to penalize gaps at the beginning of s2 seq
		matrixG[0][i] = -10000000;
	}
	double score1,score2,score3;
	time_t start = clock();
	for (int i = 1; i < matrixV.size();i++){
		for (int j = 1; j < matrixV.at(i).size(); j++){
			///V
			double prfScore = prf.getElement(i-1, s2[j][0]);
			double featPrfScore = featPrf.getScore(i-1,s2[j]);
			double gapMod = featPrf.getGapMod(i-1,s2[j]);// gap modifier, based on the features profile
			score1 = matrixV[i-1][j-1] + prfScore + featPrfScore;
			score2 = matrixG[i-1][j-1] + prfScore + featPrfScore;
			score3 = matrixH[i-1][j-1] + prfScore + featPrfScore;
			matrixV[i][j] = findVal::maxValueDoubles(score1,score2,score3);
			///G
			score1 = matrixV[i-1][j] + gapOpening * gapMod;
			score2 = matrixG[i-1][j] + gapExtension * gapMod;
			matrixG[i][j] = (score1 > score2) ? score1 : score2;
			///H
			score1 = matrixV[i][j-1] + gapOpeningHorizontal * gapMod;
			score2 = matrixH[i][j-1] + gapExtensionHorizontal * gapMod;
			matrixH[i][j] = (score1 > score2) ? score1 : score2;
		}
	}
}
//function findBestScore - returns alignment score with positions in the scoring matrix: [score, i, j] (must be either in the last row or in the last column of the scoring matrix)
std::vector<int> ScoringMatrix::findBestScore(){
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
	std::vector<int> resArr;
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
std::vector< std::vector<double> > ScoringMatrix::getVec(){
	return matrixV;
}
//function nwAlignment - performs a sequence vs profile(/pseudoprofile) needleman wunsch alignment ENCODED SEQUNCES
void ScoringMatrix::nwAlignment(std::vector<std::vector<std::string> > *result,std::vector<std::string> s2, Profile& prf, FeaturesProfile& featPrf,std::string verbose, int codon_length){
	std::vector<std::string> s1 = misc::pseudoSequence(prf.getMatrix()[0].size()+1); //creating polyA pseudoSequence representing the profile, to know later where are the gaps in the profile
	//std::string s1(prf.getMatrix()[0].size()+1,'A');
	std::string gap_code = txtProc::gapCode(codon_length);
	s2 = vecUtil::push_front(s2,gap_code);
	std::vector<std::string> newS1;
	std::vector<std::string> newS2;
	std::vector<std::vector<std::string> > ali;
	std::string newChar1;	
	std::string newChar2;	
	int i = s1.size()-1;
	int j = s2.size()-1;
	std::string currentMatrix = "V";
	int iteratorM = 10;
	//if bestScore isn't in the lower right corner, then add gaps to newS1 or newS2
	if (findBestScore().at(0) != matrixV.size()-1 || findBestScore().at(1) != matrixV.at(0).size()-1){
		i = findBestScore().at(0);
		j = findBestScore().at(1);
		for (int k = s1.size()-1; k > i; k--){
			newChar1 = s1.at(k);
			newChar2 = gap_code;
			newS1.push_back(newChar1);
			newS2.push_back(newChar2);
		}
		for (int k = s2.size()-1; k > j; k--){
			newChar1 = gap_code;
			newChar2 = s2.at(k);
			newS2.push_back(newChar2);
			newS1.push_back(newChar1);
		}
	}
	//trace back the matrix
	while (i > 0 || j > 0){
		double gapMod = featPrf.getGapMod(i-1,s2.at(j));
		if (i > 0 && j > 0 && currentMatrix == "V"){	//match/mismatch
			newChar1 = s1.at(i);
			newChar2 = s2.at(j);
			double prfScore = prf.getElement(i-1,s2.at(j)[0]);
			double featPrfScore = featPrf.getScore(i-1, s2.at(j));
			if (matrixV[i][j] != matrixV[i-1][j-1] + prfScore + featPrfScore){
				if( i > 0 && j > 0 && matrixV[i][j] == matrixG[i-1][j-1]+prfScore + featPrfScore){
					currentMatrix = "G";	
				}
				else if(i > 0 && j > 0 && matrixV[i][j] == matrixH[i-1][j-1]+prfScore + featPrfScore){
					currentMatrix = "H";
				}
			}
			i--;
			j--;
		}
		else if (i > 0 && currentMatrix == "G"){	//gap in seq2
			newChar1 = s1.at(i);
			newChar2 = gap_code;
			if (matrixG[i][j] == matrixV[i-1][j] + gapOpening * gapMod)
				currentMatrix = "V";
			i--;
		}
		else if (j > 0 && currentMatrix == "H"){	//gap in profile
			newChar1 = gap_code;
			newChar2 = s2.at(j);
			if (matrixH[i][j] == matrixV[i][j-1] + gapOpeningHorizontal * gapMod){
				currentMatrix = "V";
			}
			j--;
		}
		newS1.push_back(newChar1);
		newS2.push_back(newChar2);
	}
	if (verbose!="0"){
		misc::printEncodedSeq(newS1);	
		misc::printEncodedSeq(newS2);	
	}
	std::reverse(newS1.begin(),newS1.end());
	std::reverse(newS2.begin(),newS2.end());
	ali.push_back(newS1);
	ali.push_back(newS2);
	*result = ali;
}
