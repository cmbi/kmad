#include "misc.h"
#include "residue.h"

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>


Residue misc::CreateGapResidue(int codon_length) {
  std::string single_codon(codon_length, 'A');
  single_codon[0] = '-';
  Residue res(single_codon);
  return res;
}


ResidueSequence misc::PseudoResidueSequence(int seq_length, int codon_length) {
  std::string single_codon(codon_length, 'A');
  Residue single_res(single_codon);
  ResidueSequence result(seq_length, single_res);
  return result;
}


bool misc::CheckIfFileExists(const std::string* name) {
  const char *cstr = name->c_str();
  std::ifstream infile(cstr);
  return infile.good();
}
