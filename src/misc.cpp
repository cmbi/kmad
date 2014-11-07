#include "misc.h"
#include "Residue.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
//creates a polyA sequence of length seqLength, where each resdue is coded by a codon of length codon_length
/*
std::vector<std::string> misc::pseudoSequence(int seqLength, int codon_length){
	std::string single_codon(codon_length, 'A');
	std::vector<std::string> result(seqLength, single_codon);
	return result;
}
*/
//creates a gap codon of length codon_length
Residue misc::gapRes(int codon_length){
	std::string single_codon(codon_length,'A');
	single_codon[0] = '-';
	std::vector<std::string> additional_features;
	Residue res(single_codon, additional_features);
	return res;
}


//creates a polyA sequence of length seqLength, where each resdue is coded by a codon of length codon_length
std::vector<Residue> misc::pseudoResidueSequence(int seqLength, int codon_length){
	std::string single_codon(codon_length, 'A');
	std::vector<std::string> additional_features;	
	Residue single_res(single_codon,additional_features);
	std::vector<Residue> result(seqLength, single_res);
	return result;
}


void misc::printEncodedSeq(const std::vector<std::string>& sequence){
	for (unsigned int i = 0; i < sequence.size();i++) std::cout << sequence[i][0];	
	std::cout << "\n";
}


//look for mistakes in the given command line arguments
bool misc::checkParameters(int codonLength,int phosph,int domain,int motif,double gep,double gop,bool weightsOn, double endPenalty){
	bool alright = true;
	if (codonLength < 1 || codonLength > 10){
		alright = false;
		std::cout << "please change the codon's length to more than 0 and less than 11" << std::endl;
	}
	else if (gep >= 0 || gop >= 0 || endPenalty > 0){
		alright = false;
		std::cout << "you set gap penalty value(s) to positive (or zero)" << std::endl;
	}
	else if (phosph < 0 || domain < 0 || motif < 0){
		alright = false;
		std::cout << "you're penalizing alignment of features";
	}

	return alright;
}


bool misc::file_exists(const std::string* name){
  const char *cstr = name->c_str();
  std::ifstream infile(cstr);
  return infile.good();
}
