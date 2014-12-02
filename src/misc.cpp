#include "misc.h"
#include "residue.h"

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>


Residue misc::CreateGapResidue(int codon_length){
  std::string single_codon(codon_length,'A');
  single_codon[0] = '-';
  Residue res(single_codon);
  return res;
}


ResidueSequence misc::PseudoResidueSequence(int seq_length, int codon_length){
  std::string single_codon(codon_length, 'A');
  Residue single_res(single_codon);
  ResidueSequence result(seq_length, single_res);
  return result;
}


bool misc::CheckParameters(int codon_length, int phosph, int domain, int motif,
                           double gep, double gop, double end_pen){
  bool alright = true;
  if (codon_length < 1 || codon_length > 10){
    alright = false;
    std::cout << "please change the codon's length to more than \
                 0 and less than 11" << std::endl;
  }
  else if (gep >= 0 || gop >= 0 || end_pen > 0){
    alright = false;
    std::cout << "you set gap penalty value(s) to positive (or zero)" 
              << std::endl;
  }
  else if (phosph < 0 || domain < 0 || motif < 0){
    alright = false;
    std::cout << "you're penalizing alignment of features";
  }

  return alright;
}


bool misc::CheckIfFileExists(const std::string* name){
  const char *cstr = name->c_str();
  std::ifstream infile(cstr);
  return infile.good();
}
