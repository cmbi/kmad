#ifndef FASTA_H
#define FASTA_H

#include "types.h"

#include <iostream>
#include <fstream>
#include <vector>


namespace fasta {
  struct Residue {
    Residue(std::string codon) : codon(codon) {}

    std::string codon;
  };

  struct Sequence {
    std::string description;
    std::vector<Residue> residues;
  };

  struct FastaData {
    std::vector<fasta::Sequence> sequences;
    std::map<std::string, double> probabilities;
  };

  ///
  /// parsess a 'fasta' file, returns a Sequences object,
  /// writes motif ids and probabilities to ids and probs
  ///
  FastaData parse_fasta(std::string filename, int codonLength);

}


#endif /* FASTA_H */
