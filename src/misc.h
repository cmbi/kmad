#ifndef MISC_H
#define MISC_H

#include "types.h"
#include <iostream>
#include <vector>
#include <fstream>

class Residue;
namespace misc{
  ///
  /// creates a polyA sequence of length seqLength, where each resdue is coded 
  /// by a codon of length codon_length
  ///
  ResidueSequence PseudoResidueSequence(int seq_length, int codon_length);
}

#endif /* MISC_H */
