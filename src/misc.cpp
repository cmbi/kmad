#include "misc.h"
#include "residue.h"

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>


ResidueSequence misc::PseudoResidueSequence(int seq_length, int codon_length) {
  std::string single_codon(codon_length, 'A');
  Residue single_res(single_codon);
  ResidueSequence result(seq_length, single_res);
  return result;
}
