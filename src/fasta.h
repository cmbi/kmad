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

  typedef std::vector<fasta::Sequence> SequenceList;

  struct FastaData {
    SequenceList sequences;
    std::map<std::string, double> probabilities;
  };

  ///
  /// parsess a 'fasta' file, returns a Sequences object,
  /// writes motif ids and probabilities to ids and probs
  ///
  FastaData parse_fasta(std::string filename, int codonLength);

  Sequence make_sequence(const std::string& description,
                         const std::string& codons, int codon_length);

}


#endif /* FASTA_H */
