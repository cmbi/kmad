#include "ScoringMatrix.h"
#include "FeaturesProfile.h"
#include "Residue.h"
#include "substitutionMatrix.h"
#include "Profile.h"
#include "misc.h"
#include "vecUtil.h"
#include "txtProc.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


ScoringMatrix::ScoringMatrix(int s1size,int s2size, double pen, 
                             double endPenalty, double extensionPenalty)
:	m_iLength(s1size),
	m_jLength(s2size),	
	m_gapOpening(pen),
	m_gapExtension(extensionPenalty),
	m_endGapPenalty(endPenalty),
	m_gapOpeningHorizontal(m_gapOpening),
	m_gapExtensionHorizontal(m_gapExtension)
{
  //creates a row for the scoring matrices of length m_jLength 
  //(length of the jth sequence + 1)
	scoringMatrixRow row(m_jLength+1,0); 
  //creates a vector of vectors 'row', of length m_iLength+1 
  //(length of the ith sequence +1)
	m_matrixV.assign(m_iLength+1, row); 
	m_matrixG.assign(m_iLength+1, row);
	m_matrixH.assign(m_iLength+1, row);
}


void ScoringMatrix::calculateScores(ResidueSequence s2, Profile& prf, 
                                    FeaturesProfile& featPrf, 
                                    int codon_length){
	s2 = vecUtil::push_front(s2,misc::gapRes(codon_length));

  assert(m_matrixV.size() == m_matrixG.size());
  assert(m_matrixV.size() == m_matrixH.size());

	for (unsigned int i = 1; i < m_matrixV.size(); i++){
		m_matrixV[i][0] = -10000000; //=== infinity
		m_matrixH[i][0] = -10000000;
		m_matrixG[i][0] = i*m_endGapPenalty;	
	}
	for (unsigned int i = 1; i < m_matrixV[0].size(); i++){
		m_matrixV[0][i] = -10000000;
		m_matrixH[0][i] = i*m_endGapPenalty;	
		m_matrixG[0][i] = -10000000;
	}
	double score1,score2,score3 = 0;
	for (unsigned int i = 1; i < m_matrixV.size();i++){
		for (unsigned int j = 1; j < m_matrixV[i].size(); j++){
			///V
			double prfScore = prf.getElement(i-1, s2[j].getAA());
			double add_score = 0;
      FeaturesList features = s2[j].getFeatIndexes();
			featPrf.getScore(i-1, features, add_score);

			double final_score = prfScore + add_score;
			score1 = m_matrixV[i-1][j-1];
			score2 = m_matrixG[i-1][j-1];
			score3 = m_matrixH[i-1][j-1];

      m_matrixV[i][j] = std::max<double>(score1,
          std::max<double>(score2, score3)) + final_score;
			///G
			score1 = m_matrixV[i-1][j] + m_gapOpening;
			score2 = m_matrixG[i-1][j] + m_gapExtension;
			m_matrixG[i][j] = (score1 > score2) ? score1 : score2;
			///H
			score1 = m_matrixV[i][j-1] + m_gapOpeningHorizontal;
			score2 = m_matrixH[i][j-1] + m_gapExtensionHorizontal;
			m_matrixH[i][j] = (score1 > score2) ? score1 : score2;
		}
	}
}


valueCoords ScoringMatrix::findBestScore(){
	int maxI = m_matrixV.size()-1;
	int maxJ = m_matrixV[0].size()-1;
	int n = maxI; // last row of m_matrixV
	int m = maxJ; // last column of m_matrixV
	double maxIval = m_matrixV[maxI][maxJ];
	double maxJval = m_matrixV[maxI][maxJ];
	double real_val;
	double max = m_matrixV[maxI][maxJ];
	for (int i = 0; i < n ; i++){		//finds max score in the last row
    // add end gap penalties to the score to calc the 'real' score 
    // of the alignment
		real_val = m_matrixV[i][m]+m_endGapPenalty*(m_matrixV.size()-i); 
		if (real_val > max){
			maxIval = real_val;
			max = real_val;
			maxI = i;
		} 
	}
	for (int i = 0; i < m; i++){		//finds max score in the last column	
		real_val = m_matrixV[n][i]+m_endGapPenalty*(m_matrixV[0].size()-i); 
    // add end gap penalties to the score, to calc the 'real' score 
    // of the alignment
		if (real_val > max){
			maxJval = real_val;
			max = real_val;
			maxJ = i;
		}
	}
	valueCoords resArr;
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


void ScoringMatrix::nwAlignment(SequenceList *result,
                                ResidueSequence s2, Profile& prf, 
                                FeaturesProfile& featPrf,
                                int codon_length){
  //creating polyA pseudoSequence representing the profile, 
  //to know later where are the gaps in the profile
	ResidueSequence s1 = misc::pseudoResidueSequence(prf.getMatrix()[0].size()+1, 
                                                   codon_length); 
	Residue gap_code = misc::gapRes(codon_length);
	s2 = vecUtil::push_front(s2,gap_code);
	ResidueSequence newS1;
	ResidueSequence newS2;
	SequenceList ali; //alignment
	Residue newChar1;
	Residue newChar2;	
	int i = s1.size()-1;
	int j = s2.size()-1;
	std::string currentMatrix = "V";
	//if bestScore isn't in the lower right corner, then add gaps to newS1 or newS2
	if (findBestScore()[0] != (signed)m_matrixV.size()-1 || findBestScore()[1] != (signed)m_matrixV[0].size()-1){
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

  assert(m_matrixV.size() == m_matrixG.size());
  assert(m_matrixV.size() == m_matrixH.size());

	//trace back the matrix
	while (i > 0 || j > 0){
		//double gapMod = featPrf.getGapMod(i-1,s2[j].getFeatures());
		if (i > 0 && j > 0 && currentMatrix == "V"){	//match/mismatch
			newChar1 = s1[i];
			newChar2 = s2[j];
			double prfScore = prf.getElement(i-1,s2[j].getAA());
			double add_score = 0;
      //std::vector<std::string> features = s2[j].getFeatures();
      FeaturesList features = s2[j].getFeatIndexes();
			featPrf.getScore(i-1, features, add_score);
			double final_score = prfScore + add_score;
			if (m_matrixV[i][j] != m_matrixV[i-1][j-1] + final_score){
				if( i > 0 && j > 0 && m_matrixV[i][j] == m_matrixG[i-1][j-1]+final_score){
					currentMatrix = "G";	
				}
				else if(i > 0 && j > 0 && m_matrixV[i][j] == m_matrixH[i-1][j-1]+final_score){
					currentMatrix = "H";
				}
			}
			i--;
			j--;
		}
		else if (i > 0 && currentMatrix == "G"){	//gap in seq2
			newChar1 = s1[i];
			newChar2 = gap_code;
			if (m_matrixG[i][j] == m_matrixV[i-1][j] + m_gapOpening)
				currentMatrix = "V";
			i--;
		}
		else if (j > 0 && currentMatrix == "H"){	//gap in profile
			newChar1 = gap_code;
			newChar2 = s2[j];
			if (m_matrixH[i][j] == m_matrixV[i][j-1] + m_gapOpeningHorizontal){
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
