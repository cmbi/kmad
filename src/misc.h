#ifndef MISC_H
#define MISC_H

#include "types.h"
#include <iostream>
#include <vector>
#include <fstream>

class Residue;
namespace misc{
  ///
  /// creates a gap Residue object
  ///
	Residue CreateGapResidue(int codon_length);
  ///
  /// creates a polyA sequence of length seqLength, where each resdue is coded 
  /// by a codon of length codon_length
  ///
  ResidueSequence PseudoResidueSequence(int seq_length, int codon_length);
  ///
  /// looks for mistakes in the given command line arguments
  ///
	bool CheckParameters(int codon_length, int phosph, int domain, int motif, 
                       double gep, double gop, double end_penalty);
  ///
  /// returns True if the file 'name' exists
  ///
  bool CheckIfFileExists(const std::string* name);
}

#endif /* MISC_H */
