//UsefulStuff class implementation
#include "misc.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
//addSequenceToMultipleAlignment - throw it out?
void misc::addSequenceToMultipleAlignment( std::vector<std::string>& al1, std::vector<std::string> al2){    ///al1 - multiple alignment, al2 - pairwise alignment
	int i = 0;
	int j = 0;
	std::vector<std::string> result;
	std::string seq="";
	for (int i =0; i<al1.size()+1;i++){
		result.push_back(seq);
	}
	char gap = '-';
	while (i < al1.at(0).size() || j < al2.at(0).size()){
		char al20j = al2.at(0).at(j);
		char al10i = al1.at(0).at(i);
		if (al10i != gap && al20j != gap){
			for (int k = 0; k< al1.size();k++){
				result.at(k)+=al1.at(k).at(i);
			}
			result.at(al1.size())+=al2.at(1).at(j);		
			i++;
			j++;
		}
		else if (al10i == gap){
			for (int k = 0; k< al1.size();k++){
				result.at(k)+=al1.at(k).at(i);
			}
			result.at(al1.size())+="-";
			i++;
		}
		else{
			for (int k = 0; k< al1.size();k++){
				result.at(k)+="-";
			}
			result.at(al1.size())+=al2.at(1).at(j);
			j++;
		}
	}
	al1 = result;
}
//function countTrueValuesInVector
int misc::countTrueValuesInVector(const std::vector<bool>& vec){
	int result = 0;
	for (int i = 0; i < vec.size(); i++){
		if (vec[i]){
			result++;	
		}
	}
	return result;
}
std::vector<std::string> misc::pseudoSequence(int seqLength){
	std::vector<std::string> result(seqLength, "AAAA");
	return result;
}
void misc::printEncodedSeq(const std::vector<std::string>& sequence){
	for (int i = 0; i < sequence.size();i++) std::cout << sequence[i][0];	
	std::cout << "\n";
}
