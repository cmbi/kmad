//ScoringMatrix class implementation
//
//responsible for pairwise alignments
//creates scoring matrix, backtraces alignment path, etc.
#include "ScoringMatrix.h"
#include "FeaturesProfile.h"
#include "Residue.h"
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
ScoringMatrix::ScoringMatrix(int s1size,int s2size, double pen, 
                             double endPenalty, double extensionPenalty)
:	iLength(s1size),
	jLength(s2size),	
	gapOpening(pen),
	gapExtension(extensionPenalty),
	endGapPenalty(endPenalty),
	gapOpeningHorizontal(gapOpening),
	gapExtensionHorizontal(gapExtension)
{
	std::vector<double> row(jLength+1,0); //creates a row for the scoring matrices of length jLength (length of the jth sequence + 1)
	matrixV.assign(iLength+1, row); //creates a vector of vectors 'row', of length iLength+1 (length of the ith sequence +1)
	matrixG.assign(iLength+1, row);
	matrixH.assign(iLength+1, row);
}
//function calculateScoresProfile - calculates scoring matrix for sequences s1 and s2 using profile prf instead of a substitution matrix ENCODED SEQUENCES
void ScoringMatrix::calculateScores(std::vector<Residue> s2, Profile& prf, 
                                    FeaturesProfile& featPrf, int debug, 
                                    int codon_length, int sequence_no){
	s2 = vecUtil::push_front(s2,misc::gapRes(codon_length));
	for (unsigned int i = 1; i < matrixV.size(); i++){
		matrixV[i][0] = -10000000; //infinity
		matrixH[i][0] = -10000000;
		//matrixG[i][0] = 0;	//makes gaps at the beginning of the vertical sequence (s1) free
		matrixG[i][0] = i*endGapPenalty;	//makes gaps at the beginning of the vertical sequence (s1) free
		//matrixG[i][0] = gapOpening+(i-1)*gapExtension; // uncomment to penalize gaps at the beginnig of s1 seq with the same penalty function as the gaps inside
	}
	for (unsigned int i = 1; i < matrixV[0].size(); i++){
		matrixV[0][i] = -10000000;
		//matrixH[0][i] = 0;	//makes gaps at the beginning of the horizontal sequence (s2) free
		matrixH[0][i] = i*endGapPenalty;	//makes gaps at the beginning of the horizontal sequence (s2) free
		//matrixH[0][i] = gapOpening+(i-1)*gapExtension; //uncomment to penalize gaps at the beginning of s2 seq
		matrixG[0][i] = -10000000;
	}
	double score1,score2,score3;
	for (unsigned int i = 1; i < matrixV.size();i++){
		for (unsigned int j = 1; j < matrixV[i].size(); j++){
			///V
			double prfScore = prf.getElement(i-1, s2[j].getAA());
			double add_score = 0;
      std::vector<std::string> features = s2[j].getFeatures();
			featPrf.getScore(i-1, features, add_score);
			double final_score = prfScore + add_score;
			score1 = matrixV[i-1][j-1];
			score2 = matrixG[i-1][j-1];
			score3 = matrixH[i-1][j-1];
			matrixV[i][j] = findVal::maxValueDoubles(score1,score2,score3) + final_score;
			///G
			score1 = matrixV[i-1][j] + gapOpening;
			score2 = matrixG[i-1][j] + gapExtension;
			matrixG[i][j] = (score1 > score2) ? score1 : score2;
			///H
			score1 = matrixV[i][j-1] + gapOpeningHorizontal;
			score2 = matrixH[i][j-1] + gapExtensionHorizontal;
			matrixH[i][j] = (score1 > score2) ? score1 : score2;
		}
	}
}
//function findBestScore - returns positions of the end of the best scoring alignment[score, i, j] (must be either in the last row or in the last column of the scoring matrix)
std::vector<int> ScoringMatrix::findBestScore(){
	int maxI = matrixV.size()-1;
	int maxJ = matrixV[0].size()-1;
	int n = maxI;  //last row of matrixV
	int m = maxJ; // last column of matrixV
	double maxIval = matrixV[maxI][maxJ];
	double maxJval = matrixV[maxI][maxJ];
	double real_val;
	double max = matrixV[maxI][maxJ];
	for (int i = 0; i < n ; i++){		//finds max score in the last row
		real_val = matrixV[i][m]+endGapPenalty*(matrixV.size()-i); // add end gap penalties to the score to calc the 'real' score of the alignment
		if (real_val > max){
			maxIval = real_val;
			max = real_val;
			maxI = i;
		} 
	}
	for (int i = 0; i < m; i++){		//finds max score in the last column	
		real_val = matrixV[n][i]+endGapPenalty*(matrixV[0].size()-i); //add end gap penalties to the score, to calc the 'real' score of the alignment
		if (real_val > max){
			maxJval = real_val;
			max = real_val;
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
	return resArr; //coords of the max score
}
//function getVec - returns scoring matrix
std::vector< std::vector<double> > ScoringMatrix::getVec(){
	return matrixV;
}
//function nwAlignment - performs a sequence vs profile(/pseudoprofile) needleman wunsch alignment 
void ScoringMatrix::nwAlignment(std::vector<std::vector<Residue> > *result,
                                std::vector<Residue> s2, Profile& prf, 
                                FeaturesProfile& featPrf, std::string verbose, 
                                int codon_length, int sequence_no){
	std::vector<Residue> s1 = misc::pseudoResidueSequence(prf.getMatrix()[0].size()+1,codon_length); //creating polyA pseudoSequence representing the profile, to know later where are the gaps in the profile
	Residue gap_code = misc::gapRes(codon_length);
	s2 = vecUtil::push_front(s2,gap_code);
	std::vector<Residue> newS1;
	std::vector<Residue> newS2;
	std::vector<std::vector<Residue> > ali; //alignment
	Residue newChar1;
	Residue newChar2;	
	int i = s1.size()-1;
	int j = s2.size()-1;
	std::string currentMatrix = "V";
	//if bestScore isn't in the lower right corner, then add gaps to newS1 or newS2
	if (findBestScore()[0] != (signed)matrixV.size()-1 || findBestScore()[1] != (signed)matrixV[0].size()-1){
		i = findBestScore()[0];
		j = findBestScore()[1];
		for (int k = s1.size()-1; k > i; k--){
			newChar1 = s1[k];
			newChar2 = gap_code;
			newS1.push_back(newChar1);
			newS2.push_back(newChar2);
		}
		for (int k = s2.size()-1; k > j; k--){
			newChar1 = gap_code;
			newChar2 = s2[k];
			newS2.push_back(newChar2);
			newS1.push_back(newChar1);
		}
	}
	//trace back the matrix
	while (i > 0 || j > 0){
		//double gapMod = featPrf.getGapMod(i-1,s2[j].getFeatures());
		if (i > 0 && j > 0 && currentMatrix == "V"){	//match/mismatch
			newChar1 = s1[i];
			newChar2 = s2[j];
			double prfScore = prf.getElement(i-1,s2[j].getAA());
			double add_score = 0;
      std::vector<std::string> features = s2[j].getFeatures();
			featPrf.getScore(i-1, features, add_score);
			double final_score = prfScore + add_score;
			if (matrixV[i][j] != matrixV[i-1][j-1] + final_score){
				if( i > 0 && j > 0 && matrixV[i][j] == matrixG[i-1][j-1]+final_score){
					currentMatrix = "G";	
				}
				else if(i > 0 && j > 0 && matrixV[i][j] == matrixH[i-1][j-1]+final_score){
					currentMatrix = "H";
				}
			}
			i--;
			j--;
		}
		else if (i > 0 && currentMatrix == "G"){	//gap in seq2
			newChar1 = s1[i];
			newChar2 = gap_code;
			if (matrixG[i][j] == matrixV[i-1][j] + gapOpening)
				currentMatrix = "V";
			i--;
		}
		else if (j > 0 && currentMatrix == "H"){	//gap in profile
			newChar1 = gap_code;
			newChar2 = s2[j];
			if (matrixH[i][j] == matrixV[i][j-1] + gapOpeningHorizontal){
				currentMatrix = "V";
			}
			j--;
		}
		newS1.push_back(newChar1);
		newS2.push_back(newChar2);
	}
	std::reverse(newS1.begin(),newS1.end()); //need to reverse the sequences, because tracing back the alignment goes from the end to the beginning
	std::reverse(newS2.begin(),newS2.end());
	ali.push_back(newS1);
	ali.push_back(newS2);
	*result = ali;
}
