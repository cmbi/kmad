#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <map>
#include <vector>


namespace fasta {

  struct Residue {
    Residue(std::string codon, std::vector<std::string> features) 
      : codon(codon),
        features(features) {}
    Residue(std::string codon) : codon(codon) {}
    Residue(){};

    std::string codon;
    std::vector<std::string> features;
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
  Sequence make_sequence(const std::vector<Residue>&);
  Sequence make_sequence(unsigned long sequence_length,
                         const fasta::Residue residue);
  Sequence make_sequence(std::string sequence_string, int codon_length);
  Residue make_residue(const std::string& codon);
}


#endif /* FASTA_H */
