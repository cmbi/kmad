//UsefulStuff class implementation
#include "misc.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
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
bool misc::checkParameters(int codonLength,int phosph,int domain,int motif,double gep,double gop,bool weightsOn,double cutoff){
	bool alright = true;
	if (codonLength < 1 || codonLength > 10){
		alright = false;
		std::cout << "this is ridiculous, please change the codon's length" << std::endl;
	}
	else if (gep >= 0 || gop >= 0){
		alright = false;
		std::cout << "ohhh maaaan! you set gap penalty value(s) to positive (or zero), what were you thinking??" << std::endl;
	}
	else if (phosph < 0 || domain < 0 || motif < 0){
		alright = false;
		std::cout << "you're penalizining alignment of domains and/or phosphorylations - that's not the way to go";
	}

	return alright;
}
